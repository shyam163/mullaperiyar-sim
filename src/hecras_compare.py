"""Extract HEC-RAS 2D results and compare with our solver's matching run.

Reads hecras/run180/mulla.p01.hdf (RAS) and
outputs/baseline_142_rasmatch/ (our solver, 180 m, no lake slab), both
driven by the identical terrain + Froehlich 142 ft hydrograph.

Outputs: comparison table (stdout + hecras/comparison.md) and a
max-depth difference map (hecras/comparison_maxdepth.png).
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "hecras" / "run180"
PROJ = ROOT / "hecras" / "project"
OURS = ROOT / "outputs" / "baseline_142_rasmatch"
ARRIVE = 0.1


def load_ras():
    f = h5py.File(RUN / "mulla.p01.hdf")
    base = ("Results/Unsteady/Output/Output Blocks/Base Output/"
            "Unsteady Time Series")
    ws = f[base + "/2D Flow Areas/Perimeter 1/Water Surface"][()]
    t_days = f[base + "/Time"][()]
    t_s = (t_days - t_days[0]) * 86400.0
    g = f["Geometry/2D Flow Areas/Perimeter 1"]
    zmin = g["Cells Minimum Elevation"][()]
    cc = g["Cells Center Coordinate"][()]
    n_int = f["Geometry/2D Flow Areas/Cell Points"].shape[0]
    depth = np.clip(ws - zmin[None, :], 0, None)
    return dict(t=t_s, depth=depth[:, :n_int], cc=cc[:n_int],
                zmin=zmin[:n_int], f=f)


def series_at(ras, xy):
    i = int(np.argmin(np.linalg.norm(ras["cc"] - np.asarray(xy), axis=1)))
    return ras["depth"][:, i]


def metrics(t, d):
    if (d > ARRIVE).any():
        ta = float(t[np.argmax(d > ARRIVE)])
        return ta / 3600, float(d.max()), float(t[d.argmax()]) / 3600
    return None, 0.0, None


def main():
    ras = load_ras()
    gauges = json.loads((PROJ / "gauges.json").read_text())

    ours_meta = json.load(open(OURS / "run_meta.json"))
    ours = np.load(OURS / "results.npz")
    gs = ours["gauge_series"]
    names = ours_meta["gauge_names"]
    t_ours = np.arange(gs.shape[0]) * 60.0

    rows = []
    for g in ["dam_toe", "vandiperiyar"]:
        d_ras = series_at(ras, gauges[g])
        a1, p1, tp1 = metrics(ras["t"], d_ras)
        d_our = gs[:, names.index(g)]
        a2, p2, tp2 = metrics(t_ours, d_our)
        rows.append((g, a2, p2, a1, p1))

    # basin (Idukki) volume gain over time in both models
    import rasterio
    with rasterio.open(ROOT / "data" / "dem_utm180.tif") as src:
        tr = src.transform
    br, bc_ = gauges["_basin_rows_cols"]
    bxy = np.c_[tr.c + (np.array(bc_) + 0.5) * tr.a,
                tr.f + (np.array(br) + 0.5) * tr.e]
    # RAS basin cells: cells within 200 m of any (subsampled) basin cell
    from scipy.spatial import cKDTree
    tree = cKDTree(bxy)
    dist, _ = tree.query(ras["cc"], k=1)
    bsel = dist < 300.0
    area_per_cell = 180.0 * 180.0
    vol_ras = (ras["depth"][:, bsel] * area_per_cell).sum(axis=1)
    a_ras = None
    if (vol_ras > 1e4).any():
        a_ras = float(ras["t"][np.argmax(vol_ras > 1e4)]) / 3600
    a_our = (ours_meta["sink_first_arrival_s"] or 0) / 3600

    # mass check for RAS: final volume vs hydrograph total
    t5, q5 = np.loadtxt(PROJ / "hydrograph_5min.csv", delimiter=",",
                        skiprows=1).T
    v_in = np.trapezoid(q5, t5)
    v_ras_end = float((ras["depth"][-1] * area_per_cell).sum())
    print(f"RAS volume: in {v_in/1e6:.0f} Mm3, stored at end "
          f"{v_ras_end/1e6:.0f} Mm3 ({v_ras_end/v_in*100:.1f}%)")

    lines = [
        "# HEC-RAS 2D vs our solver - same terrain, same hydrograph",
        "",
        "180 m mesh/grid, Froehlich 142 ft hydrograph, Idukki dry-start,",
        "spillway closed. HEC-RAS 6.6 (full SWE, implicit, PARDISO);",
        "ours = explicit Rusanov FV (numba).",
        "",
        "| gauge | ours: arrival (h) | ours: peak (m) | RAS: arrival (h) | RAS: peak (m) |",
        "|---|---|---|---|---|",
    ]
    for g, a2, p2, a1, p1 in rows:
        lines.append(f"| {g} | {a2:.2f} | {p2:.1f} | "
                     f"{a1 if a1 is None else round(a1,2)} | {p1:.1f} |")
    lines.append(f"| idukki basin (>1e4 m3) | {a_our:.2f} | - | "
                 f"{a_ras if a_ras is None else round(a_ras,2)} | - |")
    lines += ["",
              f"- RAS stored volume at t=24 h: {v_ras_end/1e6:.0f} of "
              f"{v_in/1e6:.0f} Mm3 injected ({v_ras_end/v_in*100:.1f}%)",
              f"- our ledger error: "
              f"{ours_meta['ledger_error']*100:+.3f}%"]
    out = "\n".join(lines)
    print(out)
    (ROOT / "hecras" / "comparison.md").write_text(out + "\n")

    # ---- max depth difference map on our grid
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    md_ras_cells = ras["depth"].max(axis=0)
    with rasterio.open(ROOT / "data" / "dem_utm180.tif") as src:
        shape = (src.height, src.width)
        tr = src.transform
    cols = np.clip(((ras["cc"][:, 0] - tr.c) / tr.a - 0.5).round(), 0,
                   shape[1] - 1).astype(int)
    rows_ = np.clip(((ras["cc"][:, 1] - tr.f) / tr.e - 0.5).round(), 0,
                    shape[0] - 1).astype(int)
    ras_md = np.zeros(shape, np.float32)
    ras_md[rows_, cols] = md_ras_cells
    our_md = ours["max_depth"]
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    for ax, dat, ttl in [(axes[0], our_md, "ours: max depth"),
                         (axes[1], ras_md, "HEC-RAS 2D: max depth"),
                         (axes[2], ras_md - our_md, "RAS - ours")]:
        cmap, vmin, vmax = ("Blues", 0, 20) if "max" in ttl else \
            ("RdBu", -8, 8)
        im = ax.imshow(np.where(np.abs(dat) > 0.05, dat, np.nan),
                       cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, shrink=0.7)
        ax.set_title(ttl)
        ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(ROOT / "hecras" / "comparison_maxdepth.png", dpi=80)
    print("wrote hecras/comparison.md + comparison_maxdepth.png")


if __name__ == "__main__":
    main()
