"""Compare LISFLOOD-FP results with our solver's matching run.

Both engines: identical 180 m conditioned terrain (dam wall + Idukki
blocks), identical Froehlich 142 ft hydrograph at the same gorge cells,
identical Manning map, dry Idukki start, 24 h.
LISFLOOD-FP = local-inertial (acceleration) solver;
ours = explicit full-SWE Rusanov FV.

Writes lisflood/comparison.md and lisflood/comparison_maxdepth.png.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
LF = ROOT / "lisflood"
OURS = ROOT / "outputs" / "baseline_142_rasmatch"
ARRIVE = 0.1


def read_asc(path):
    with open(path) as f:
        hdr = {}
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
        arr = np.loadtxt(f)
    return arr, hdr


def metrics(t, d):
    if (d > ARRIVE).any():
        return (float(t[np.argmax(d > ARRIVE)]) / 3600, float(np.max(d)),
                float(t[np.argmax(d)]) / 3600)
    return None, 0.0, None


def main():
    # ---- lisflood outputs
    names = json.loads((LF / "gauge_names.json").read_text())
    stage = np.loadtxt(LF / "results" / "mulla.stage", skiprows=10)
    t_lf = stage[:, 0]
    # .stage reports water surface elevation? no: stage output is depth?
    # LISFLOOD stage file columns: time, then stage (water depth) per gauge
    lf_series = {n: stage[:, i + 1] for i, n in enumerate(names)}
    md_lf, hdr = read_asc(LF / "results" / "mulla.max")
    md_lf[md_lf < 0] = 0.0

    mass = np.loadtxt(LF / "results" / "mulla.mass", skiprows=1)
    t_m, vol, qin = mass[:, 0], mass[:, 5], mass[:, 6]
    v_in_total = np.trapezoid(qin, t_m)

    # ---- our matching run
    meta = json.load(open(OURS / "run_meta.json"))
    res = np.load(OURS / "results.npz")
    gs = res["gauge_series"]
    gnames = meta["gauge_names"]
    t_our = np.arange(gs.shape[0]) * 60.0
    md_our = res["max_depth"]

    # basin volume arrival for lisflood: use idukki_pool gauge + basin mask
    # from our domain (same grid) applied to hourly depth rasters is
    # overkill; the pool gauge arrival is the comparable signal.
    rows = []
    for g in ["dam_toe", "vandiperiyar", "idukki_pool"]:
        a_lf, p_lf, tp_lf = metrics(t_lf, lf_series[g])
        d_our = gs[:, gnames.index(g)]
        a_our, p_our, tp_our = metrics(t_our, d_our)
        rows.append((g, a_our, p_our, tp_our, a_lf, p_lf, tp_lf))

    fmt = lambda v: "-" if v is None else f"{v:.2f}"
    lines = [
        "# LISFLOOD-FP vs our solver - same terrain, same hydrograph",
        "",
        "180 m grid, Froehlich 142 ft hydrograph at the same gorge cells,",
        "dry Idukki start, spillway closed, 24 h. LISFLOOD-FP 8.0.3",
        "local-inertial ('acceleration') vs our explicit full-SWE Rusanov FV.",
        "",
        "| gauge | ours: arrival (h) | ours: peak (m) | ours: t_peak (h) |"
        " LF: arrival (h) | LF: peak (m) | LF: t_peak (h) |",
        "|---|---|---|---|---|---|---|",
    ]
    for g, a2, p2, tp2, a1, p1, tp1 in rows:
        lines.append(f"| {g} | {fmt(a2)} | {p2:.1f} | {fmt(tp2)} | "
                     f"{fmt(a1)} | {p1:.1f} | {fmt(tp1)} |")
    wet_lf = int((md_lf > ARRIVE).sum())
    wet_our = int((md_our > ARRIVE).sum())
    both = ((md_lf > ARRIVE) & (md_our > ARRIVE)).sum()
    iou = both / max(((md_lf > ARRIVE) | (md_our > ARRIVE)).sum(), 1)
    lines += [
        "",
        f"- inflow volume (LF mass file): {v_in_total/1e6:.0f} Mm3 "
        f"(hydrograph 380)",
        f"- final stored volume (LF): {vol[-1]/1e6:.0f} Mm3",
        f"- wet footprint >10 cm: ours {wet_our} cells, LF {wet_lf} cells, "
        f"IoU {iou:.2f}",
        f"- our ledger error {meta['ledger_error']*100:+.3f}%; "
        f"LF cumulative Verror {np.abs(mass[:,10]).sum()/1e6:.3f} Mm3",
    ]
    out = "\n".join(lines)
    print(out)
    (LF / "comparison.md").write_text(out + "\n")

    # ---- max depth maps
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    for ax, dat, ttl, cmap, vmin, vmax in [
            (axes[0], md_our, "ours (Rusanov FV): max depth", "Blues", 0, 20),
            (axes[1], md_lf, "LISFLOOD-FP: max depth", "Blues", 0, 20),
            (axes[2], md_lf - md_our, "LF - ours [m]", "RdBu", -8, 8)]:
        im = ax.imshow(np.where(np.abs(dat) > 0.05, dat, np.nan),
                       cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, shrink=0.7)
        ax.set_title(ttl)
        ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(LF / "comparison_maxdepth.png", dpi=80)
    print("wrote lisflood/comparison.md + comparison_maxdepth.png")


if __name__ == "__main__":
    main()
