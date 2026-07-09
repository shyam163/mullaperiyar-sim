"""Author a headless HEC-RAS 2D project for the Mullaperiyar baseline.

Generates everything the hecras-v66-linux toolchain's Workflow C needs:
    hecras/project/
        mulla.prj  mulla.g01  mulla.p01  mulla.u01
        Terrain/terrain.tif
        Land Classification/LandCover.{tif,hdf}
        gauges.json  hydrograph_5min.csv

Same physics inputs as our solver: the terrain with dam wall + Idukki
blocks burned in, the Froehlich 142 ft hydrograph on a BC line across the
gorge below the dam, Manning 0.06/0.035 split at 300 m elevation.
The 2D mesh covers only the flood corridor (baseline max-depth footprint
dilated ~2 km), not the whole domain.

Usage: python src/export_hecras.py [--res 180]
"""
import argparse
import json
import sys
import uuid
from pathlib import Path

import h5py
import numpy as np
import rasterio
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import terrain  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT / "hecras" / "project"
NAME = "mulla"
AREA = "Perimeter 1"      # ras_preprocess hardcodes this in Plan Parameters


def f16(v):
    """One coordinate as an exactly-16-char field (HEC fixed format)."""
    s = f"{v:.7f}"
    if len(s) > 16:
        s = s[:16]
    return s.rjust(16)


def corridor_polygon(dom, dx):
    """Flood-corridor perimeter from the baseline run footprint."""
    from shapely.geometry import Polygon
    from skimage import measure

    res = np.load(ROOT / "outputs" / "baseline_142" / "results.npz")
    md = res["max_depth"]
    # our baseline ran at 90 m; resample its footprint onto this grid if
    # needed (nearest is fine for a mask)
    if md.shape != dom.z.shape:
        zoom = dom.z.shape[0] / md.shape[0]
        md = ndimage.zoom(md, zoom, order=0)
        md = md[:dom.z.shape[0], :dom.z.shape[1]]
    mask = md > 0.05
    mask[dom.sink_rows, dom.sink_cols] = True
    # the reservoir is NOT part of the RAS domain (the hydrograph replaces
    # it); drop its footprint so the corridor perimeter crosses the gorge
    # at the dam and the inflow BC line can sit on the perimeter there
    resmask = np.zeros_like(mask)
    resmask[dom.res_rows, dom.res_cols] = True
    mask &= ~ndimage.binary_dilation(resmask, iterations=3)
    lab, _ = ndimage.label(ndimage.binary_dilation(mask, iterations=3))
    keep = lab[dom.inj_rows[0], dom.inj_cols[0]]
    mask = lab == keep
    mask = ndimage.binary_dilation(mask, iterations=int(round(2000.0 / dx)))
    mask = ndimage.binary_fill_holes(mask)
    # upstream cut: remove everything on the reservoir side of the
    # injection point within 6 km of the dam, so the perimeter crosses
    # the gorge right at the injection cells
    dam_rc = dom.ll_to_rc(*terrain.DAM_MULLA)
    inj_r, inj_c = float(dom.inj_rows.mean()), float(dom.inj_cols.mean())
    u = np.array([dam_rc[0] - inj_r, dam_rc[1] - inj_c])
    u /= np.linalg.norm(u) + 1e-9
    rr, cc = np.mgrid[:mask.shape[0], :mask.shape[1]]
    upstream = ((rr - inj_r) * u[0] + (cc - inj_c) * u[1]) > 1.0
    near = ((rr - dam_rc[0]) ** 2 + (cc - dam_rc[1]) ** 2) < (
        6000.0 / dx) ** 2
    mask &= ~(upstream & near)
    lab, _ = ndimage.label(mask)
    mask = lab == lab[int(round(inj_r - 3 * u[0])), int(round(inj_c - 3 * u[1]))]

    cont = max(measure.find_contours(mask.astype(float), 0.5), key=len)
    # rows/cols -> UTM
    tr = dom.transform
    xs = tr.c + (cont[:, 1] + 0.5) * tr.a
    ys = tr.f + (cont[:, 0] + 0.5) * tr.e
    poly = Polygon(np.c_[xs, ys]).simplify(1.5 * dx)
    poly = poly.buffer(0)
    x, y = poly.exterior.xy
    return np.c_[x, y], mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, default=180)
    args = ap.parse_args()
    S = float(args.res)

    (PROJ / "Terrain").mkdir(parents=True, exist_ok=True)
    (PROJ / "Land Classification").mkdir(parents=True, exist_ok=True)

    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(ROOT / "data" / f"dem_utm{args.res}.tif",
                                   380e6)
    tr = dom.transform
    wkt = rasterio.crs.CRS.from_epsg(32643).to_wkt()

    # ---------------- terrain (walls + blocks burned in dom.z)
    # RAS subgrid hydraulic tables need terrain much FINER than the mesh
    # (one pixel per cell -> degenerate volume curves -> solver segfault),
    # so export at 30 m: a bilinear upsample of the conditioned 90 m grid.
    from rasterio.warp import Resampling, reproject
    with rasterio.open(ROOT / "data" / "dem_utm90.tif") as src90:
        import contextlib as _ctx
        import io as _io
        with _ctx.redirect_stdout(_io.StringIO()):
            dom90 = terrain.build_domain(ROOT / "data" / "dem_utm90.tif",
                                         380e6)
        tr90 = src90.transform
    tr30 = rasterio.Affine(30.0, 0, tr90.c, 0, -30.0, tr90.f)
    h30, w30 = dom90.z.shape[0] * 3, dom90.z.shape[1] * 3
    z30 = np.empty((h30, w30), np.float32)
    reproject(dom90.z.astype(np.float32), z30, src_transform=tr90,
              src_crs="EPSG:32643", dst_transform=tr30, dst_crs="EPSG:32643",
              resampling=Resampling.bilinear)
    prof = dict(driver="GTiff", height=h30, width=w30, count=1,
                dtype="float32", crs="EPSG:32643", transform=tr30,
                compress="deflate")
    with rasterio.open(PROJ / "Terrain" / "terrain.tif", "w", **prof) as d:
        d.write(z30, 1)

    # ---------------- land cover: 1 = gorge forest, 2 = plains
    lc = np.where(dom.z > 300.0, 1, 2).astype(np.uint8)
    prof["dtype"] = "uint8"
    with rasterio.open(PROJ / "Land Classification" / "LandCover.tif", "w",
                       **prof) as d:
        d.write(lc, 1)
    with h5py.File(PROJ / "Land Classification" / "LandCover.hdf", "w") as f:
        f.attrs["File Type"] = np.bytes_("HEC Land Cover")
        f.attrs["GUID"] = np.bytes_(str(uuid.uuid4()))
        f.attrs["LC Type"] = np.bytes_("LandCover")
        f.attrs["Projection"] = np.bytes_(wkt)
        f.attrs["Version"] = np.bytes_("2.0")
        rm = np.array([(0, b"NoData"), (1, b"ForestGorge"), (2, b"Plains")],
                      dtype=[("ID", "<i4"), ("Name", "S28")])
        f.create_dataset("Raster Map", data=rm)
        va = np.array([(b"NoData", 0.052, 0.0),
                       (b"ForestGorge", 0.0892, 0.0),
                       (b"Plains", 0.052, 0.0)],
                      dtype=[("Name", "S28"), ("ManningsN", "<f4"),
                             ("Percent Impervious", "<f4")])
        f.create_dataset("Variables", data=va)

    # ---------------- corridor mesh: perimeter + seed points
    perim, cmask = corridor_polygon(dom, S)
    from matplotlib.path import Path as MplPath
    ny, nx = dom.z.shape
    xs = tr.c + (np.arange(nx) + 0.5) * tr.a
    ys = tr.f + (np.arange(ny) + 0.5) * tr.e
    step = int(round(S / abs(tr.a)))
    gx, gy = np.meshgrid(xs[::step], ys[::step])
    pts = np.c_[gx.ravel(), gy.ravel()]
    inside = MplPath(perim).contains_points(pts, radius=-0.3 * S)
    seeds = pts[inside]
    print(f"corridor: {len(perim)} perimeter pts, {len(seeds)} seed cells "
          f"@ {S:.0f} m")

    # ---------------- BC lines
    # inflow: across the gorge at the injection cells, perpendicular-ish;
    # just use a segment centered on the mean injection point, oriented
    # along the line joining the two outermost injection cells rotated 90deg
    inj_xy = np.c_[tr.c + (dom.inj_cols + 0.5) * tr.a,
                   tr.f + (dom.inj_rows + 0.5) * tr.e]
    c0 = inj_xy.mean(axis=0)
    d = inj_xy[-1] - inj_xy[0]
    d = d / (np.linalg.norm(d) + 1e-9)
    n = np.array([-d[1], d[0]])
    L = max(3.5 * S, 400.0)
    bc_in = (c0 - n * L, c0 + n * L)
    # outlet (inactive, stays dry): far southwest corner of the corridor,
    # snapped onto high ground inside the perimeter
    sw = perim[np.argmin(perim[:, 0] + perim[:, 1])]
    bc_out = (sw + np.array([S, 2 * S]), sw + np.array([4 * S, 2 * S]))

    # ---------------- .g01
    g = []
    g.append(f"Geom Title={NAME}")
    g.append("Program Version=6.60")
    g.append(f"Storage Area={AREA:<16s},{c0[0]:.7f},{c0[1]:.7f}")
    g.append(f"Storage Area Surface Line= {len(perim)} ")
    for x, y in perim:
        g.append(f16(x) + f16(y) + " " * 16)
    g.append("Storage Area Type= 1 ")
    g.append("Storage Area Area=")
    g.append("Storage Area Min Elev=")
    g.append("Storage Area Is2D=-1")
    g.append(f"Storage Area Point Generation Data=,,{S:.0f},{S:.0f}")
    g.append(f"Storage Area 2D Points= {len(seeds)} ")
    for i in range(0, len(seeds), 2):
        row = "".join(f16(v) for v in seeds[i:i + 2].ravel())
        g.append(row)
    g.append("Storage Area 2D PointsPerimeterTime=01Jan2026 00:00:00")
    g.append("Storage Area Mannings=0.0892")
    for nm, (a, b) in [("Inflow", bc_in), ("Outlet", bc_out)]:
        mid = (np.asarray(a) + np.asarray(b)) / 2
        g.append(f"BC Line Name={nm:<32s}")
        g.append(f"BC Line Storage Area={AREA:<16s}")
        g.append(f"BC Line Start Position= {a[0]:.7f} , {a[1]:.7f} ")
        g.append(f"BC Line Middle Position= {mid[0]:.7f} , {mid[1]:.7f} ")
        g.append(f"BC Line End Position= {b[0]:.7f} , {b[1]:.7f} ")
        g.append("BC Line Arc= 2 ")
        g.append(f16(a[0]) + f16(a[1]) + f16(b[0]) + f16(b[1]))
        g.append(f"BC Line Text Position= {mid[0]:.7f} , {mid[1]:.7f} ")
    (PROJ / f"{NAME}.g01").write_text("\n".join(g) + "\n")

    # ---------------- hydrograph at 5-min resolution [m3/s]
    t, q = np.load(ROOT / "outputs" / "baseline_142" / "hydrograph.npy")
    t5 = np.arange(0, 86400 + 1, 300.0)
    q5 = np.interp(t5, t, q)
    np.savetxt(PROJ / "hydrograph_5min.csv",
               np.c_[t5, q5], delimiter=",", header="t_s,q_m3s")

    # ---------------- .u01
    u = []
    u.append(f"Flow Title={NAME}")
    u.append("Program Version=6.60")
    u.append("Use Restart= 0 ")
    u.append(f"Boundary Location=                ,                ,        "
             f",        ,                ,{AREA:<16s},                "
             f",Inflow                          ,                ")
    u.append("Interval=5MIN")
    u.append(f"Flow Hydrograph= {len(q5)} ")
    for i in range(0, len(q5), 10):
        u.append("".join(f"{v:8.1f}" for v in q5[i:i + 10]))
    u.append("Stage Hydrograph TW Check=0")
    u.append("Flow Hydrograph Slope= 0.005 ")
    u.append("Use DSS=False")
    u.append("Use Fixed Start Time=False")
    u.append("Fixed Start Date/Time=,")
    u.append("Is Critical Boundary=False")
    u.append("Critical Boundary Flow=")
    u.append(f"Boundary Location=                ,                ,        "
             f",        ,                ,{AREA:<16s},                "
             f",Outlet                          ,                ")
    u.append("Friction Slope=0.001,0")
    (PROJ / f"{NAME}.u01").write_text("\n".join(u) + "\n")

    # ---------------- .p01
    p = f"""Plan Title={NAME}
Program Version=6.60
Short Identifier=01
Simulation Date=01JAN2026,0000,02JAN2026,0000
Geom File=g01
Flow File=u01
Computation Interval=10SEC
Output Interval=5MIN
Instantaneous Interval=5MIN
Mapping Interval=5MIN
Run HTab=-1
Run UNet=-1
Run PostProcess= 0
UNET Theta=1
UNET ZTol=0.003
UNET ZSATol=0.003
UNET MxIter=20
UNET DZMax Abort=30
UNET D2 Equation=1
UNET D2 Theta=1
UNET D2 Theta Warmup=1
UNET D2 Z Tol=0.003
UNET D2 Volume Tol=0.003
UNET D2 Max Iterations=20
UNET D2 TimeSlices=1
UNET D2 Cores=16
UNET D2 Turbulence Formulation=None
UNET D2 SolverType=PARDISO (Direct)
Write IC File= 0
UNET Use Existing IC File= 0
"""
    (PROJ / f"{NAME}.p01").write_text(p)

    # ---------------- .prj
    (PROJ / f"{NAME}.prj").write_text(
        f"Proj Title={NAME}\nCurrent Plan=p01\nDefault Exp/Contr=0.3,0.1\n"
        f"SI Units\nGeom File=g01\nUnsteady File=u01\nPlan File=p01\n"
        f"Y Axis Title=Elevation\nX Axis Title(PF)=Main Channel Distance\n")

    # ---------------- gauges for extraction
    gauges = {}
    marks = dict(terrain.TOWNS)
    for nm in ["vandiperiyar"]:
        r, c = dom.gauges[nm]
        gauges[nm] = [tr.c + (c + 0.5) * tr.a, tr.f + (r + 0.5) * tr.e]
    r, c = int(dom.inj_rows[0]), int(dom.inj_cols[0])
    gauges["dam_toe"] = [tr.c + (c + 0.5) * tr.a, tr.f + (r + 0.5) * tr.e]
    rc = dom.ll_to_rc(9.820, 76.940)
    gauges["idukki_pool"] = [tr.c + (rc[1] + 0.5) * tr.a,
                             tr.f + (rc[0] + 0.5) * tr.e]
    # basin cell list (for volume-arrival metric), as UTM bounding info
    gauges["_basin_rows_cols"] = [dom.sink_rows.tolist()[::20],
                                  dom.sink_cols.tolist()[::20]]
    gauges["_bc_inflow"] = [list(map(float, bc_in[0])),
                            list(map(float, bc_in[1]))]
    gauges["_bc_outlet"] = [list(map(float, bc_out[0])),
                            list(map(float, bc_out[1]))]
    (PROJ / "gauges.json").write_text(json.dumps(gauges))
    print("project written to", PROJ)


if __name__ == "__main__":
    main()
