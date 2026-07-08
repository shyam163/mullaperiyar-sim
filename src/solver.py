"""2D shallow-water solver: first-order finite volume on a regular grid.

Numerics
--------
* Rusanov (local Lax-Friedrichs) flux with Audusse et al. (2004) hydrostatic
  reconstruction, which makes the scheme well-balanced: a lake at rest over
  arbitrary bathymetry is an exact steady state.
* Adaptive timestep at CFL 0.4 on the fastest wave (|u| + sqrt(g h)).
* Wet/dry handling with a thin-film threshold H_EPS: cells below it hold
  mass but carry no momentum, and faces between two dry cells are skipped.
* Manning friction as a semi-implicit (unconditionally stable) source term.
* Two numba passes per step - flux then update - each parallel over rows
  and writing to disjoint arrays, so there are no write races.

The Simulation driver adds breach-hydrograph injection, reservoir drawdown,
an "infinite sink" mask (Idukki reservoir), open outflow edges (west/south),
gauges, snapshots, and a running mass ledger.
"""
import json
import time
from pathlib import Path

import numpy as np

try:
    from numba import njit, prange
    NUMBA = True
except ImportError:  # pragma: no cover - spec requires a numpy fallback
    NUMBA = False
    prange = range

    def njit(*args, **kwargs):
        def wrap(f):
            return f
        return wrap if not (args and callable(args[0])) else args[0]

G = 9.81
H_EPS = 0.01      # thin-film threshold [m]: below this, no momentum
V_CAP = 35.0      # safety cap on speed [m/s] (gorge flows peak ~20-30)
ARRIVAL_DEPTH = 0.1  # depth defining "flood arrival" [m]


# ----------------------------------------------------------------------
# kernels
# ----------------------------------------------------------------------
@njit(parallel=True, fastmath=True, cache=True)
def kernel_dt(h, hu, hv, i0, i1, j0, j1, dx, cfl, dt_max):
    ny = h.shape[0]
    row_min = np.full(ny, dt_max)
    for i in prange(i0, i1 + 1):
        m = dt_max
        for j in range(j0, j1 + 1):
            hij = h[i, j]
            if hij > H_EPS:
                u = hu[i, j] / hij
                v = hv[i, j] / hij
                c = np.sqrt(G * hij)
                s = max(abs(u), abs(v)) + c
                if s > 1e-8:
                    d = cfl * dx / s
                    if d < m:
                        m = d
        row_min[i] = m
    return row_min[i0:i1 + 1].min()


@njit(inline="always")
def _face_flux(hL, uL, vL, zL, hR, uR, vR, zR, out):
    """Rusanov flux across one face (normal = 'x' by convention; the caller
    swaps u/v for y-faces). out[0:3] = (F_h, F_hun, F_hut) where 'un' is the
    normal momentum and 'ut' the transverse; out[3], out[4] = reconstructed
    depths on the L and R side (needed for the well-balanced bed source)."""
    zf = max(zL, zR)
    hLs = max(hL + zL - zf, 0.0)
    hRs = max(hR + zR - zf, 0.0)
    out[3] = hLs
    out[4] = hRs
    if hLs <= 0.0 and hRs <= 0.0:
        out[0] = 0.0
        out[1] = 0.0
        out[2] = 0.0
        return
    cL = np.sqrt(G * hLs)
    cR = np.sqrt(G * hRs)
    a = max(abs(uL) + cL, abs(uR) + cR)
    fL0 = hLs * uL
    fL1 = hLs * uL * uL + 0.5 * G * hLs * hLs
    fL2 = hLs * uL * vL
    fR0 = hRs * uR
    fR1 = hRs * uR * uR + 0.5 * G * hRs * hRs
    fR2 = hRs * uR * vR
    out[0] = 0.5 * (fL0 + fR0) - 0.5 * a * (hRs - hLs)
    out[1] = 0.5 * (fL1 + fR1) - 0.5 * a * (hRs * uR - hLs * uL)
    out[2] = 0.5 * (fL2 + fR2) - 0.5 * a * (hRs * vR - hLs * vL)


@njit(parallel=True, fastmath=True, cache=True)
def kernel_fluxes(z, h, hu, hv, Fx, Fy, i0, i1, j0, j1):
    """Fill face fluxes for the active window.

    Fx[i, j] = flux across the face between cells (i, j-1) and (i, j),
    for j in [j0, j1+1]. Fy[i, j] = flux between (i-1, j) and (i, j),
    for i in [i0, i1+1]. Requires 1 <= i0, i1 <= ny-2 (same for j).
    """
    # x-direction faces
    for i in prange(i0, i1 + 1):
        for j in range(j0, j1 + 2):
            hL = h[i, j - 1]
            hR = h[i, j]
            if hL <= H_EPS and hR <= H_EPS:
                # still must zero (stale values from a previous window)
                Fx[i, j, 0] = 0.0
                Fx[i, j, 1] = 0.0
                Fx[i, j, 2] = 0.0
                Fx[i, j, 3] = 0.0
                Fx[i, j, 4] = 0.0
                continue
            uL = hu[i, j - 1] / hL if hL > H_EPS else 0.0
            vL = hv[i, j - 1] / hL if hL > H_EPS else 0.0
            uR = hu[i, j] / hR if hR > H_EPS else 0.0
            vR = hv[i, j] / hR if hR > H_EPS else 0.0
            _face_flux(hL, uL, vL, z[i, j - 1], hR, uR, vR, z[i, j],
                       Fx[i, j])
    # y-direction faces (normal velocity is v)
    for i in prange(i0, i1 + 2):
        for j in range(j0, j1 + 1):
            hL = h[i - 1, j]
            hR = h[i, j]
            if hL <= H_EPS and hR <= H_EPS:
                Fy[i, j, 0] = 0.0
                Fy[i, j, 1] = 0.0
                Fy[i, j, 2] = 0.0
                Fy[i, j, 3] = 0.0
                Fy[i, j, 4] = 0.0
                continue
            uL = hu[i - 1, j] / hL if hL > H_EPS else 0.0
            vL = hv[i - 1, j] / hL if hL > H_EPS else 0.0
            uR = hu[i, j] / hR if hR > H_EPS else 0.0
            vR = hv[i, j] / hR if hR > H_EPS else 0.0
            # swap: normal = v, transverse = u
            _face_flux(hL, vL, uL, z[i - 1, j], hR, vR, uR, z[i, j],
                       Fy[i, j])


@njit(parallel=True, fastmath=True, cache=True)
def kernel_update(z, h, hu, hv, Fx, Fy, gn2, max_depth, arrival, t_now,
                  i0, i1, j0, j1, dx, dt):
    """Apply flux divergence + well-balanced bed source, then friction and
    wet/dry cleanup. Also maintains max-depth and arrival-time rasters.
    Returns the (tiny) mass clipped from negative depths, for the ledger."""
    rdx = dt / dx
    clipped = 0.0
    for i in prange(i0, i1 + 1):
        local_clip = 0.0
        for j in range(j0, j1 + 1):
            # skip cells with no water and no flux on any face
            fw = Fx[i, j, 0]
            fe = Fx[i, j + 1, 0]
            fn = Fy[i, j, 0]
            fs = Fy[i + 1, j, 0]
            hij = h[i, j]
            if (hij <= 0.0 and fw == 0.0 and fe == 0.0 and fn == 0.0
                    and fs == 0.0):
                continue

            hn = hij - rdx * (fe - fw) - rdx * (fs - fn)
            # momentum: flux divergence + Audusse interface source terms
            # (the g/2 h^2 pieces reconstruct the bed-slope force exactly
            # at rest; see module docstring)
            hun = (hu[i, j]
                   - rdx * (Fx[i, j + 1, 1] - Fx[i, j, 1])
                   - rdx * (Fy[i + 1, j, 2] - Fy[i, j, 2])
                   + rdx * 0.5 * G * (Fx[i, j + 1, 3] ** 2
                                      - Fx[i, j, 4] ** 2))
            hvn = (hv[i, j]
                   - rdx * (Fx[i, j + 1, 2] - Fx[i, j, 2])
                   - rdx * (Fy[i + 1, j, 1] - Fy[i, j, 1])
                   + rdx * 0.5 * G * (Fy[i + 1, j, 3] ** 2
                                      - Fy[i, j, 4] ** 2))

            if hn < 0.0:
                local_clip += -hn
                hn = 0.0
            if hn <= H_EPS:
                hun = 0.0
                hvn = 0.0
            else:
                # semi-implicit Manning friction
                u = hun / hn
                v = hvn / hn
                sp = np.sqrt(u * u + v * v)
                denom = 1.0 + dt * gn2[i, j] * sp / hn ** (4.0 / 3.0)
                hun /= denom
                hvn /= denom
                # safety cap
                u = hun / hn
                v = hvn / hn
                sp = np.sqrt(u * u + v * v)
                if sp > V_CAP:
                    scale = V_CAP / sp
                    hun *= scale
                    hvn *= scale

            h[i, j] = hn
            hu[i, j] = hun
            hv[i, j] = hvn
            if hn > max_depth[i, j]:
                max_depth[i, j] = hn
            if hn > ARRIVAL_DEPTH and arrival[i, j] < 0.0:
                arrival[i, j] = t_now
        clipped += local_clip
    return clipped * dx * dx


@njit(parallel=True, cache=True)
def kernel_wet_bbox(h, i0, i1, j0, j1):
    """Bounding box of wet cells inside [i0-.., expanded search window]."""
    ny, nx = h.shape
    lo = max(i0, 0)
    hi = min(i1, ny - 1)
    row_has = np.zeros(ny, np.int64)
    row_jmin = np.full(ny, nx, np.int64)
    row_jmax = np.full(ny, -1, np.int64)
    for i in prange(lo, hi + 1):
        jmn, jmx = nx, -1
        for j in range(max(j0, 0), min(j1, nx - 1) + 1):
            if h[i, j] > 0.0:
                if j < jmn:
                    jmn = j
                if j > jmx:
                    jmx = j
        if jmx >= 0:
            row_has[i] = 1
            row_jmin[i] = jmn
            row_jmax[i] = jmx
    imn, imx, jmn, jmx = ny, -1, nx, -1
    for i in range(lo, hi + 1):
        if row_has[i]:
            if i < imn:
                imn = i
            if i > imx:
                imx = i
            if row_jmin[i] < jmn:
                jmn = row_jmin[i]
            if row_jmax[i] > jmx:
                jmx = row_jmax[i]
    return imn, imx, jmn, jmx


@njit(parallel=True, cache=True)
def kernel_sum(h, i0, i1, j0, j1):
    tot = 0.0
    for i in prange(i0, i1 + 1):
        s = 0.0
        for j in range(j0, j1 + 1):
            s += h[i, j]
        tot += s
    return tot


# ----------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------
class Simulation:
    """Owns the state arrays and runs the time loop for one scenario."""

    def __init__(self, z, dx, n_map, out_dir,
                 duration=86400.0, snap_dt=300.0, cfl=0.4, dt_max=10.0):
        ny, nx = z.shape
        self.z = np.ascontiguousarray(z, np.float64)
        self.dx = float(dx)
        self.cell_area = self.dx * self.dx
        self.gn2 = np.ascontiguousarray(G * n_map ** 2, np.float64)
        self.h = np.zeros((ny, nx))
        self.hu = np.zeros((ny, nx))
        self.hv = np.zeros((ny, nx))
        self.Fx = np.zeros((ny, nx + 1, 5))
        self.Fy = np.zeros((ny + 1, nx, 5))
        self.max_depth = np.zeros((ny, nx), np.float32)
        self.arrival = np.full((ny, nx), -1.0, np.float32)
        self.duration = duration
        self.snap_dt = snap_dt
        self.cfl = cfl
        self.dt_max = dt_max
        self.out_dir = Path(out_dir)
        (self.out_dir / "snapshots").mkdir(parents=True, exist_ok=True)

        # sources/sinks configured by the scenario
        self.injections = []   # list of dicts: rows, cols, t, q, drain cells
        self.gauges = {}       # name -> (row, col)
        # monitored basin (Idukki): the pool rises physically behind the
        # sealed rim (spillway assumed CLOSED - no releases); we track its
        # net volume gain and fire the cascade hook when the surge arrives
        self.basin_cells = None
        self.on_basin_arrival = None
        self.basin_trigger_vol = 1e6
        self._basin_cb_fired = False

        # mass ledger [m3]
        self.led = dict(initial=0.0, injected=0.0, drained=0.0, sink=0.0,
                        boundary=0.0, clipped=0.0)

    # -------------------------------------------------- configuration
    def set_initial_water(self, depth):
        self.h[:] = depth
        self.led["initial"] = float(depth.sum()) * self.cell_area

    def add_injection(self, rows, cols, t, q, drain_cells=None):
        """Inject hydrograph q(t) at cells; optionally draw the same volume
        down from `drain_cells` (the upstream pool) for visual/mass
        consistency - the hydrograph remains authoritative if they run dry.
        """
        self.injections.append(dict(
            rows=np.asarray(rows), cols=np.asarray(cols),
            t=np.asarray(t, float), q=np.asarray(q, float),
            drain=(np.asarray(drain_cells[0]), np.asarray(drain_cells[1]))
            if drain_cells is not None else None))

    def set_basin(self, rows, cols):
        self.basin_cells = (np.asarray(rows), np.asarray(cols))

    # -------------------------------------------------- helpers
    def _window_from_wet(self):
        ny, nx = self.h.shape
        i0, i1, j0, j1 = kernel_wet_bbox(self.h, 0, ny - 1, 0, nx - 1)
        if i1 < 0:  # nothing wet yet
            i0, i1, j0, j1 = 1, 2, 1, 2
        for inj in self.injections:
            i0 = min(i0, inj["rows"].min())
            i1 = max(i1, inj["rows"].max())
            j0 = min(j0, inj["cols"].min())
            j1 = max(j1, inj["cols"].max())
        m = 4  # margin so the window never trails the front
        return (max(i0 - m, 1), min(i1 + m, ny - 2),
                max(j0 - m, 1), min(j1 + m, nx - 2))

    def storage(self):
        return float(self.h.sum()) * self.cell_area

    def ledger_error(self):
        led = self.led
        # clipping negative depths to zero ADDS mass, so it enters with +
        expected = (led["initial"] + led["injected"] - led["drained"]
                    - led["sink"] - led["boundary"] + led["clipped"])
        denom = max(led["initial"] + led["injected"], 1.0)
        return (self.storage() - expected) / denom

    # -------------------------------------------------- main loop
    def run(self, log_every=3600.0, gauge_dt=60.0):
        t = 0.0
        step = 0
        next_snap = 0.0
        next_log = 0.0
        n_bins = int(self.duration / gauge_dt) + 2
        gauge_names = list(self.gauges)
        gauge_rows = np.array([self.gauges[g][0] for g in gauge_names], int)
        gauge_cols = np.array([self.gauges[g][1] for g in gauge_names], int)
        gauge_series = np.zeros((n_bins, len(gauge_names)), np.float32)
        sink_bins = np.zeros(n_bins, np.float64)  # m3 per gauge_dt bin
        snap_times = []
        sink_hit_t = None
        # basin gain is measured against the initial pool (nominal slab)
        basin_vol0 = 0.0
        if self.basin_cells is not None:
            rr, cc = self.basin_cells
            basin_vol0 = float(self.h[rr, cc].sum()) * self.cell_area
        basin_prev = basin_vol0
        basin_peak = 0.0
        wall0 = time.time()

        win = self._window_from_wet()
        while t < self.duration:
            i0, i1, j0, j1 = win
            dt = kernel_dt(self.h, self.hu, self.hv, i0, i1, j0, j1,
                           self.dx, self.cfl, self.dt_max)
            dt = min(dt, self.duration - t + 1e-9)

            kernel_fluxes(self.z, self.h, self.hu, self.hv,
                          self.Fx, self.Fy, i0, i1, j0, j1)
            clip = kernel_update(self.z, self.h, self.hu, self.hv,
                                 self.Fx, self.Fy, self.gn2,
                                 self.max_depth, self.arrival, t,
                                 i0, i1, j0, j1, self.dx, dt)
            self.led["clipped"] += clip

            # breach inflow (and matched upstream-pool drawdown)
            for inj in self.injections:
                q = np.interp(t, inj["t"], inj["q"], left=0.0, right=0.0)
                if q <= 0.0:
                    continue
                vol = q * dt
                self.h[inj["rows"], inj["cols"]] += vol / (
                    self.cell_area * len(inj["rows"]))
                self.led["injected"] += vol
                if inj["drain"] is not None:
                    rr, cc = inj["drain"]
                    hr = self.h[rr, cc]
                    wet = hr > 0.0
                    nwet = int(wet.sum())
                    if nwet:
                        dh = vol / (self.cell_area * nwet)
                        removed = np.minimum(hr, dh * wet)
                        self.h[rr, cc] = hr - removed
                        self.led["drained"] += float(removed.sum()) * \
                            self.cell_area

            # Idukki basin monitor: the pool rises physically (rim sealed,
            # spillway closed); track net volume gain for reporting and
            # fire the cascade hook when the surge has clearly arrived
            if self.basin_cells is not None:
                rr, cc = self.basin_cells
                vol = float(self.h[rr, cc].sum()) * self.cell_area
                gain_step = vol - basin_prev
                basin_prev = vol
                if gain_step > 0.0:
                    b = min(int(t / gauge_dt), n_bins - 1)
                    sink_bins[b] += gain_step
                gain = vol - basin_vol0
                basin_peak = max(basin_peak, gain)
                if sink_hit_t is None and gain >= 1e4:
                    sink_hit_t = t
                if (self.on_basin_arrival and not self._basin_cb_fired
                        and gain >= self.basin_trigger_vol):
                    self._basin_cb_fired = True
                    self.on_basin_arrival(t)

            # open boundaries: 2-cell strips on all four edges swallow
            # outflow (kernel never updates the outermost cells, so keeping
            # these strips empty guarantees zero flux across the domain rim
            # and makes every exit show up in the ledger)
            for sl in (np.s_[:, 0], np.s_[:, 1], np.s_[:, -1], np.s_[:, -2],
                       np.s_[0, :], np.s_[1, :], np.s_[-1, :], np.s_[-2, :]):
                hb = self.h[sl]
                vol = float(hb.sum()) * self.cell_area
                if vol > 0.0:
                    self.led["boundary"] += vol
                    self.h[sl] = 0.0
                    self.hu[sl] = 0.0
                    self.hv[sl] = 0.0

            # gauges: keep the max depth seen inside each 60 s bin
            b = min(int(t / gauge_dt), n_bins - 1)
            np.maximum(gauge_series[b], self.h[gauge_rows, gauge_cols],
                       out=gauge_series[b])

            t += dt
            step += 1

            if t >= next_snap:
                idx = len(snap_times)
                np.savez_compressed(
                    self.out_dir / "snapshots" / f"h_{idx:04d}.npz",
                    h_cm=np.clip(self.h * 100.0, 0, 65535).astype(np.uint16),
                    t=t)
                snap_times.append(t)
                next_snap += self.snap_dt

            # active window: CFL < 1 means the front moves at most one cell
            # per step, so growing by 1 cell keeps the window airtight;
            # a periodic full scan re-tightens it.
            if step % 25 == 0:
                win = self._window_from_wet()
            else:
                ny, nx = self.h.shape
                win = (max(win[0] - 1, 1), min(win[1] + 1, ny - 2),
                       max(win[2] - 1, 1), min(win[3] + 1, nx - 2))

            if t >= next_log:
                err = self.ledger_error()
                print(f"t={t/3600:7.3f} h  step={step:8d}  dt={dt:6.3f} s  "
                      f"wet={(self.h > H_EPS).sum():7d}  "
                      f"stor={self.storage()/1e6:8.1f} Mm3  "
                      f"ledger_err={err*100:+.3f}%  "
                      f"wall={time.time()-wall0:7.1f} s", flush=True)
                next_log += log_every

        # ---------------- final bookkeeping
        err = self.ledger_error()
        result = dict(
            duration=self.duration, dx=self.dx, steps=step,
            wall_s=time.time() - wall0,
            ledger={k: v for k, v in self.led.items()},
            ledger_error=err,
            storage_final=self.storage(),
            sink_total=basin_peak,   # peak volume impounded in the basin
            sink_first_arrival_s=sink_hit_t,
            snap_times=snap_times,
            gauge_names=gauge_names,
            gauge_dt=gauge_dt,
        )
        np.savez_compressed(
            self.out_dir / "results.npz",
            max_depth=self.max_depth, arrival=self.arrival,
            gauge_series=gauge_series, sink_bins=sink_bins,
            snap_times=np.array(snap_times))
        with open(self.out_dir / "run_meta.json", "w") as f:
            json.dump(result, f, indent=2, default=float)
        print(f"DONE steps={step} wall={result['wall_s']:.0f}s "
              f"ledger_err={err*100:+.3f}%"
              + ("  [FLAG >2%]" if abs(err) > 0.02 else ""), flush=True)
        return result
