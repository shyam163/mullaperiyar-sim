"""Domain setup on the real DEM: masks, walls, sources, gauges.

Everything geographic lives here so the solver can stay abstract.
Grid convention: row 0 = north (rasterio north-up), col 0 = west.
"""
from dataclasses import dataclass

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform
from scipy import ndimage

# --- locations (lat, lon) -------------------------------------------------
DAM_MULLA = (9.5286, 77.1394)
DAM_CHERUTHONI = (9.8450, 76.9770)
TOWNS = {
    "vandiperiyar": (9.566, 77.088),
    "neriamangalam": (10.05, 76.78),
    "kalady": (10.17, 76.44),
    "aluva": (10.11, 76.35),
    "varappuzha": (10.07, 76.27),
}
DAM_IDUKKI_ARCH = (9.8425, 76.9760)
DAM_KULAMAVU = (9.8117, 76.8936)
RES_BOX = (9.44, 9.58, 77.08, 77.26)       # Mullaperiyar reservoir clip box
RES_SEED = (9.531, 77.170)                 # inside Thekkady lake
IDUKKI_BOX = (9.70, 9.95, 76.85, 77.13)    # Idukki reservoir clip box
# several candidate seeds spread over the lake; only those whose DSM value
# looks like water surface (690-736 m) are used, and the fills are unioned.
# (a single snap-to-lowest seed can land in the gorge downstream instead)
IDUKKI_SEEDS = [(9.838, 76.980), (9.820, 76.940), (9.780, 77.040),
                (9.800, 76.990), (9.760, 77.080), (9.830, 77.060)]
IDUKKI_MAX_WSE = 736.0                     # FRL 732.6 m + margin
IDUKKI_SURFACE = (690.0, 736.0)            # plausible DSM lake-surface range
import os
LAKE_SLAB = float(os.environ.get("MULLA_LAKE_SLAB", 5.0))  # initial lake depth [m]


@dataclass
class Domain:
    z: np.ndarray            # bed elevation [m], possibly with dam wall burned
    dx: float
    transform: rasterio.Affine
    crs: object
    n_map: np.ndarray        # Manning coefficient
    res_rows: np.ndarray     # burned reservoir cells
    res_cols: np.ndarray
    depth0: np.ndarray       # initial water depth
    wse: float               # solved reservoir surface elevation
    inj_rows: np.ndarray     # Mullaperiyar breach injection cells
    inj_cols: np.ndarray
    sink_rows: np.ndarray    # Idukki reservoir sink cells
    sink_cols: np.ndarray
    cheru_rows: np.ndarray   # Cheruthoni breach injection cells
    cheru_cols: np.ndarray
    gauges: dict             # name -> (row, col)

    def ll_to_rc(self, lat, lon):
        xs, ys = warp_transform("EPSG:4326", self.crs, [lon], [lat])
        r, c = rasterio.transform.rowcol(self.transform, xs[0], ys[0])
        return int(r), int(c)

    def rc_to_ll(self, r, c):
        x, y = rasterio.transform.xy(self.transform, r, c)
        lon, lat = warp_transform(self.crs, "EPSG:4326", [x], [y])
        return lat[0], lon[0]


def _box_slices(dom_transform, crs, box, shape):
    """(lat0, lat1, lon0, lon1) -> row/col slices on the UTM grid."""
    lat0, lat1, lon0, lon1 = box
    xs, ys = warp_transform("EPSG:4326", crs,
                            [lon0, lon1, lon0, lon1],
                            [lat0, lat0, lat1, lat1])
    rows, cols = [], []
    for x, y in zip(xs, ys):
        r, c = rasterio.transform.rowcol(dom_transform, x, y)
        rows.append(r)
        cols.append(c)
    r0, r1 = max(min(rows), 0), min(max(rows), shape[0] - 1)
    c0, c1 = max(min(cols), 0), min(max(cols), shape[1] - 1)
    return slice(r0, r1 + 1), slice(c0, c1 + 1)


def _connected_below(z, wse, seed_rc, rslice, cslice):
    """Cells with z < wse, 4-connected to the seed, inside the clip box."""
    sub = z[rslice, cslice] < wse
    lab, _ = ndimage.label(sub)
    sr, sc = seed_rc[0] - rslice.start, seed_rc[1] - cslice.start
    if not (0 <= sr < sub.shape[0] and 0 <= sc < sub.shape[1]) or not sub[sr, sc]:
        return None
    comp = lab == lab[sr, sc]
    rows, cols = np.nonzero(comp)
    return rows + rslice.start, cols + cslice.start


def _annulus_lowest(z, center_rc, dx, r_in, r_out, bearing_deg, half_angle,
                    k, exclude=None):
    """K lowest-elevation cells in an annulus, restricted to a bearing cone.

    bearing: 0 = north, 90 = east (grid: row decreases northwards).
    """
    r0, c0 = center_rc
    nrad = int(np.ceil(r_out / dx)) + 1
    cand = []
    for dr in range(-nrad, nrad + 1):
        for dc in range(-nrad, nrad + 1):
            r, c = r0 + dr, c0 + dc
            if not (0 <= r < z.shape[0] and 0 <= c < z.shape[1]):
                continue
            dist = np.hypot(dr, dc) * dx
            if not (r_in <= dist <= r_out):
                continue
            ang = (np.degrees(np.arctan2(dc, -dr)) + 360.0) % 360.0
            dang = min(abs(ang - bearing_deg), 360 - abs(ang - bearing_deg))
            if dang > half_angle:
                continue
            if exclude is not None and exclude[r, c]:
                continue
            cand.append((z[r, c], r, c))
    cand.sort()
    picked = cand[:k]
    # drop cells far above the channel floor (hillside cells would just
    # dribble the hydrograph down a slope)
    zmin = picked[0][0]
    picked = [p for p in picked if p[0] <= zmin + 25.0]
    return (np.array([p[1] for p in picked]),
            np.array([p[2] for p in picked]))


def build_domain(dem_path, target_volume_m3):
    """Load the DEM and assemble everything a scenario needs."""
    with rasterio.open(dem_path) as src:
        z = src.read(1).astype(np.float64)
        transform_ = src.transform
        crs = src.crs
        dx = abs(transform_.a)

    dom = Domain(z=z, dx=dx, transform=transform_, crs=crs,
                 n_map=None, res_rows=None, res_cols=None, depth0=None,
                 wse=0.0, inj_rows=None, inj_cols=None, sink_rows=None,
                 sink_cols=None, cheru_rows=None, cheru_cols=None, gauges={})

    def plug(latlon, radius_m, crest_m):
        """Raise a disc of cells to the dam crest: the 30 m DSM smooths
        narrow concrete crests below pool level, so flood fills would
        otherwise leak straight through the dam."""
        rc = dom.ll_to_rc(*latlon)
        rad = int(np.ceil(radius_m / dx))
        rr, cc = np.ogrid[:z.shape[0], :z.shape[1]]
        disc = ((rr - rc[0]) ** 2 + (cc - rc[1]) ** 2) <= rad ** 2
        z[disc] = np.maximum(z[disc], crest_m)
        return disc

    def plug_rect(box, crest_m):
        """Solid rectangular wall block. Disc plugs and rim rings kept
        leaking through the DSM's smoothed crests/spillway notches; a
        solid block spanning the whole dam site cannot have gaps."""
        rsl_, csl_ = _box_slices(transform_, crs, box, z.shape)
        blk = np.zeros(z.shape, bool)
        blk[rsl_, csl_] = z[rsl_, csl_] < crest_m
        z[rsl_, csl_] = np.maximum(z[rsl_, csl_], crest_m)
        return blk

    # ---- dam walls
    dam_rc = dom.ll_to_rc(*DAM_MULLA)
    wall = plug(DAM_MULLA, 400.0, 930.0)
    # Idukki arch + Cheruthoni + the NE corridor exit, as one block along
    # the dam line; Kulamavu saddle as another. Crest 742 m (~dam crests).
    idukki_plugs = (plug_rect((9.8405, 9.8555, 76.9410, 77.0280), 742.0)
                    | plug_rect((9.7960, 9.8125, 76.8760, 76.9000), 742.0))

    # ---- reservoir: bisect WSE until the connected fill matches volume
    rsl, csl = _box_slices(transform_, crs, RES_BOX, z.shape)
    seed = dom.ll_to_rc(*RES_SEED)
    # make sure the seed is on the lake bed (lowest DSM cell nearby)
    sr = slice(max(seed[0] - 5, 0), seed[0] + 6)
    sc = slice(max(seed[1] - 5, 0), seed[1] + 6)
    off = np.unravel_index(np.argmin(z[sr, sc]), z[sr, sc].shape)
    seed = (sr.start + off[0], sc.start + off[1])

    lake_floor = z[seed]
    lo, hi = lake_floor + 0.5, lake_floor + 80.0
    cells = None
    for _ in range(48):
        wse = 0.5 * (lo + hi)
        got = _connected_below(z, wse, seed, rsl, csl)
        if got is None:
            lo = wse
            continue
        rows, cols = got
        vol = float(np.sum(wse - z[rows, cols])) * dx * dx
        if vol < target_volume_m3:
            lo = wse
        else:
            hi = wse
            cells = (rows, cols, wse, vol)
    if cells is None:
        raise RuntimeError("reservoir fill failed")
    rows, cols, wse, vol = cells
    depth0 = np.zeros_like(z)
    depth0[rows, cols] = wse - z[rows, cols]
    dom.res_rows, dom.res_cols, dom.wse = rows, cols, wse
    dom.depth0 = depth0
    area_km2 = len(rows) * dx * dx / 1e6
    print(f"reservoir: WSE={wse:.1f} m ASL, {vol/1e6:.0f} Mm3, "
          f"{area_km2:.1f} km2, floor={lake_floor:.0f} m")

    # ---- Mullaperiyar breach injection cells: gorge floor NW of the dam
    resmask = np.zeros(z.shape, bool)
    resmask[rows, cols] = True
    resmask |= wall
    dom.inj_rows, dom.inj_cols = _annulus_lowest(
        z, dam_rc, dx, 350.0, 900.0, bearing_deg=315.0, half_angle=90.0,
        k=6, exclude=resmask)
    print(f"injection cells at z={z[dom.inj_rows, dom.inj_cols].mean():.0f} m,"
          f" {len(dom.inj_rows)} cells")

    # ---- Idukki sink mask: union of fills from verified lake-surface seeds
    isl, icsl = _box_slices(transform_, crs, IDUKKI_BOX, z.shape)
    sinkmask = np.zeros(z.shape, bool)
    used = 0
    for lat, lon in IDUKKI_SEEDS:
        iseed = dom.ll_to_rc(lat, lon)
        zval = z[iseed]
        if not (IDUKKI_SURFACE[0] <= zval <= IDUKKI_SURFACE[1]):
            continue
        got = _connected_below(z, IDUKKI_MAX_WSE, iseed, isl, icsl)
        if got is None:
            continue
        sinkmask[got[0], got[1]] = True
        used += 1
    if used == 0:
        raise RuntimeError("no valid Idukki lake seed found")
    sinkmask &= ~idukki_plugs
    # The 30 m DSM smooths the dam crests and the Kulamavu saddle below
    # pool level, so the fill leaks into the downstream gorge no matter
    # how carefully the dams are plugged. The lake surface itself sits at
    # 724-733 m and entirely south of the dam line, while the leak paths
    # drop below ~695 m or run north of the dams - so keep only
    # lake-plausible cells:
    sinkmask &= (z >= 695.0) & (z <= IDUKKI_MAX_WSE)
    r_northlimit, _ = dom.ll_to_rc(9.848, 76.977)
    sinkmask[:r_northlimit, :] = False

    # ---- verify the basin is watertight: fill from the lake at just
    # below block crest; the fill must stay inside the reservoir's
    # geographic footprint (it may run up the inflow arms - legitimate
    # backwater - but never north past the dams or west past Kulamavu)
    wide = _box_slices(transform_, crs, (9.60, 10.00, 76.80, 77.25),
                       z.shape)
    vseed = None
    for lat, lon in IDUKKI_SEEDS:
        rc_ = dom.ll_to_rc(lat, lon)
        if IDUKKI_SURFACE[0] <= z[rc_] <= IDUKKI_SURFACE[1]:
            vseed = rc_
            break
    got = _connected_below(z, 741.5, vseed, *wide)
    la0, lo0 = dom.rc_to_ll(int(got[0].min()), int(got[1].min()))
    la1, lo1 = dom.rc_to_ll(int(got[0].max()), int(got[1].max()))
    print(f"idukki basin verify: fill extent lat {la1:.3f}..{la0:.3f} "
          f"lon {lo0:.3f}..{lo1:.3f}")
    if la0 > 9.87 or lo0 < 76.85:
        raise RuntimeError("Idukki basin leaks past the dam blocks")
    dom.sink_rows, dom.sink_cols = np.nonzero(sinkmask)
    area = len(dom.sink_rows) * dx * dx / 1e6
    print(f"idukki sink: {used} seeds, {len(dom.sink_rows)} cells, "
          f"{area:.0f} km2")
    # nominal FLAT pool on the lake (surface = 725 m datum + LAKE_SLAB)
    # so an arriving surge levels out at gravity-wave speed instead of
    # crawling over a dry flat; the real bathymetry below the DSM water
    # surface is unknowable anyway. A flat surface is hydrostatically
    # consistent (well-balanced scheme keeps it exactly at rest); a
    # uniform-depth blanket over the sloping shoreline would avalanche.
    pool_wse = 725.0 + LAKE_SLAB
    lake_depth = np.maximum(pool_wse - z[dom.sink_rows, dom.sink_cols], 0.0)
    depth0[dom.sink_rows, dom.sink_cols] = lake_depth
    print(f"idukki initial pool: surface {pool_wse:.0f} m ASL, "
          f"{float(lake_depth.sum())*dx*dx/1e6:.0f} Mm3 nominal slab")
    if not (15.0 <= area <= 120.0):
        print(f"WARNING: Idukki sink area {area:.0f} km2 implausible "
              f"(FRL area ~60 km2)")

    # ---- Cheruthoni cascade injection: the gorge runs NORTH of the dams
    # (the lake sits on the south side; verified against the DSM). The
    # annulus starts beyond the rim shell so the hydrograph lands in the
    # open gorge, not on the sealed rim.
    cheru_rc = dom.ll_to_rc(*DAM_CHERUTHONI)
    dom.cheru_rows, dom.cheru_cols = _annulus_lowest(
        z, cheru_rc, dx, 1400.0, 2800.0,
        bearing_deg=0.0, half_angle=85.0,
        k=6, exclude=sinkmask | idukki_plugs)
    print(f"cheruthoni injection at "
          f"z={z[dom.cheru_rows, dom.cheru_cols].mean():.0f} m")

    # ---- Manning roughness: forested gorge above 300 m, plains below
    # experiment knob: carve a synthetic channel (one cell wide ~= 100 m,
    # MULLA_CARVE metres deep) along the steepest-descent path from the
    # dam to the Idukki basin, compensating the canopy DSM's canyon loss.
    carve_depth = float(os.environ.get("MULLA_CARVE", 0.0))
    carved = np.zeros(z.shape, bool)
    if carve_depth > 0.0:
        r, c = int(dom.inj_rows[0]), int(dom.inj_cols[0])
        visited = set()
        path = []
        for _ in range(5000):
            path.append((r, c))
            visited.add((r, c))
            if sinkmask[r, c] or z[r, c] < 700.0:
                break
            best, bestz = None, 1e18
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    rr, cc2 = r + dr, c + dc
                    if (dr == 0 and dc == 0) or (rr, cc2) in visited:
                        continue
                    if not (0 <= rr < z.shape[0] and 0 <= cc2 < z.shape[1]):
                        continue
                    if z[rr, cc2] < bestz:
                        bestz, best = z[rr, cc2], (rr, cc2)
            if best is None:
                break
            r, c = best
        cut = [(r, c) for (r, c) in path if not sinkmask[r, c]]
        mode = os.environ.get("MULLA_CARVE_MODE", "grade")
        if mode == "sill":
            # sill removal: bed = running minimum of the DSM along the
            # thalweg. Cuts only through the spurious canopy sills,
            # preserves every natural steep drop, adds minimal storage.
            bed = z[cut[0][0], cut[0][1]]
            for (r, c) in cut:
                bed = min(bed, z[r, c])
                z[r, c] = max(bed, 726.0)
                carved[r, c] = True
            print(f"synthetic channel (sill removal): {len(cut)} cells, "
                  f"end bed {z[cut[-1]]:.0f} m")
        else:
            # graded channel: linear profile from (start - carve_depth)
            # down to the lake entry (726 m) - ~ the river's mean gradient
            d = np.zeros(len(cut))
            for k in range(1, len(cut)):
                d[k] = d[k - 1] + dx * np.hypot(
                    cut[k][0] - cut[k - 1][0], cut[k][1] - cut[k - 1][1])
            z0 = z[cut[0]] - carve_depth
            z1 = 726.0
            grade = z0 + (z1 - z0) * d / max(d[-1], 1.0)
            for k, (r, c) in enumerate(cut):
                z[r, c] = min(z[r, c], grade[k])   # never raise terrain
                carved[r, c] = True
            print(f"synthetic channel: {len(cut)} cells, graded "
                  f"{z0:.0f} -> {z1:.0f} m over {d[-1]/1000:.1f} km "
                  f"(slope {(z0-z1)/max(d[-1],1):.4f})")

    n_gorge = float(os.environ.get("MULLA_N_GORGE", 0.06))
    dom.n_map = np.where(z > 300.0, n_gorge, 0.035)
    dom.n_map[carved] = 0.035   # rock/water channel, not forest

    # ---- gauges: snap each town to the lowest cell within 300 m
    snap = int(np.ceil(300.0 / dx))
    for name, (lat, lon) in TOWNS.items():
        r, c = dom.ll_to_rc(lat, lon)
        sr = slice(max(r - snap, 0), r + snap + 1)
        sc = slice(max(c - snap, 0), c + snap + 1)
        off = np.unravel_index(np.argmin(z[sr, sc]), z[sr, sc].shape)
        dom.gauges[name] = (sr.start + off[0], sc.start + off[1])

    # sanity: nothing wet may sit on the outflow rim (2-cell edge strips)
    ny, nx = z.shape
    on_rim = ((rows <= 1) | (rows >= ny - 2) | (cols <= 1) | (cols >= nx - 2))
    assert not on_rim.any(), "reservoir touches domain rim"
    return dom
