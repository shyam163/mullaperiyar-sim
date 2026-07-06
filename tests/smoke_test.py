"""Solver smoke tests on a tiny synthetic valley (no DEM needed).

1. Lake at rest over bumpy bathymetry: velocities must stay ~0 and the
   surface flat (well-balancedness of the hydrostatic reconstruction).
2. Closed-basin dam break: mass conserved to round-off, no negative
   depths, and the wet front actually advances across the dry bed.
3. Injection + ledger: hydrograph volume in == stored volume (closed box).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from solver import Simulation, H_EPS  # noqa: E402

OUT = Path(__file__).resolve().parent / "_smoke_out"


def make_sim(z, dx=10.0, n=0.0, duration=600.0):
    n_map = np.full(z.shape, n)
    sim = Simulation(z, dx, n_map, OUT, duration=duration, snap_dt=1e9,
                     dt_max=1.0)
    return sim


def test_lake_at_rest():
    ny, nx = 60, 80
    rng = np.random.default_rng(42)
    # bumpy closed basin: high walls, random interior bathymetry
    z = rng.uniform(0.0, 4.0, (ny, nx))
    z[[0, 1, -2, -1], :] = 50.0
    z[:, [0, 1, -2, -1]] = 50.0
    wse = 5.0
    sim = make_sim(z, duration=300.0)
    depth = np.maximum(wse - z, 0.0)
    depth[z >= 50.0] = 0.0
    sim.set_initial_water(depth)
    v0 = sim.storage()
    sim.run(log_every=1e9)
    umax = np.abs(sim.hu).max() + np.abs(sim.hv).max()
    surf = np.where(sim.h > 0, sim.h + sim.z, np.nan)
    flat = np.nanmax(surf) - np.nanmin(surf)
    dm = abs(sim.storage() - v0) / v0
    print(f"lake-at-rest: max|hu| {umax:.2e}, surface range {flat:.2e} m, "
          f"mass drift {dm:.2e}")
    assert umax < 1e-8, "well-balancedness violated"
    assert dm < 1e-12, "mass not conserved at rest"
    return True


def test_dam_break():
    ny, nx = 60, 200
    x = np.arange(nx) * 10.0
    # V-shaped valley sloping down +x, high rim so the basin is closed
    yy = np.abs(np.arange(ny) - ny / 2) * 0.5
    z = yy[:, None] + np.maximum(0.0, 20.0 - 0.01 * x)[None, :]
    z[[0, 1, -2, -1], :] = 100.0
    z[:, [0, 1, -2, -1]] = 100.0
    sim = make_sim(z, duration=600.0)
    # column of water behind a "dam" at x index 40
    depth = np.zeros((ny, nx))
    behind = np.s_[2:-2, 2:40]
    depth[behind] = np.maximum(30.0 - z[behind], 0.0)
    sim.set_initial_water(depth)
    v0 = sim.storage()
    front0 = np.max(np.nonzero(sim.h.max(axis=0) > H_EPS)[0])
    sim.run(log_every=1e9)
    front1 = np.max(np.nonzero(sim.h.max(axis=0) > H_EPS)[0])
    dm = abs(sim.storage() + sim.led["boundary"] + sim.led["clipped"] - v0) / v0
    neg = float(sim.h.min())
    print(f"dam-break: front {front0} -> {front1}, mass err {dm:.2e}, "
          f"min h {neg:.2e}, clipped {sim.led['clipped']:.3e} m3")
    assert front1 > front0 + 50, "wet front did not advance"
    assert dm < 1e-10, "mass not conserved in dam break"
    assert neg >= 0.0, "negative depth"
    return True


def test_injection_ledger():
    ny, nx = 60, 80
    z = np.zeros((ny, nx))
    z[[0, 1, -2, -1], :] = 50.0
    z[:, [0, 1, -2, -1]] = 50.0
    sim = make_sim(z, duration=600.0, n=0.03)
    t = np.array([0.0, 100.0, 300.0, 1e9])
    q = np.array([0.0, 500.0, 0.0, 0.0])
    sim.add_injection([30, 30], [40, 41], t, q)
    sim.run(log_every=1e9)
    err = abs(sim.ledger_error())
    print(f"injection: injected {sim.led['injected']:.1f} m3, "
          f"stored {sim.storage():.1f} m3, ledger err {err:.2e}")
    assert err < 1e-10, "ledger inconsistent"
    return True


if __name__ == "__main__":
    ok = all([test_lake_at_rest(), test_dam_break(), test_injection_ledger()])
    print("ALL SMOKE TESTS PASSED" if ok else "FAILURES")
