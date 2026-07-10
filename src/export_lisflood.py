"""Author a LISFLOOD-FP 8 project for the Mullaperiyar baseline.

Same physics inputs as our solver's matching run (180 m, no Idukki slab):
identical conditioned terrain with dam wall + Idukki blocks, identical
Froehlich 142 ft hydrograph injected at the gorge cells below the dam,
identical Manning map. Solver: LISFLOOD-FP's local-inertial
("acceleration") scheme.

Point-source units note (iterateq.cpp:500): QVAR values are m^2/s
(discharge per metre width); the model applies H += q*dx*dt/dA, i.e.
Q = q*dx. We therefore write q = Q/(dx*n_points) at each of the
injection cells.

Writes hecras/../lisflood/ project dir. Usage:
    MULLA_LAKE_SLAB=0 python src/export_lisflood.py
"""
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import terrain  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "lisflood"
RES = 180


def write_asc(path, arr, tr, fmt="%.3f"):
    ny, nx = arr.shape
    hdr = (f"ncols {nx}\nnrows {ny}\nxllcorner {tr.c:.3f}\n"
           f"yllcorner {tr.f + tr.e * ny:.3f}\ncellsize {tr.a:.3f}\n"
           f"NODATA_value -9999\n")
    with open(path, "w") as f:
        f.write(hdr)
        np.savetxt(f, arr, fmt=fmt)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", choices=["baseline", "cascade"],
                    default="baseline")
    args = ap.parse_args()
    cascade = args.scenario == "cascade"
    out = OUT / "cascade" if cascade else OUT
    out.mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(exist_ok=True)

    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(ROOT / "data" / f"dem_utm{RES}.tif", 380e6)
    tr = dom.transform
    if cascade:
        # match our cascade_142_lfmatch run: slab pool + prescribed trigger
        assert terrain.LAKE_SLAB > 0.0, "cascade needs the default slab"
        meta = json.load(open(ROOT / "outputs" / "cascade_142_lfmatch" /
                              "scenario_meta.json"))
        trigger_t = float(meta["cascade"]["trigger_t"])
        print(f"prescribing Cheruthoni breach at t={trigger_t/3600:.2f} h "
              f"(from our 180 m run)")
    else:
        assert terrain.LAKE_SLAB == 0.0, "baseline runs with MULLA_LAKE_SLAB=0"

    write_asc(out / "dem.asc", dom.z, tr)
    write_asc(out / "manning.asc", dom.n_map, tr, fmt="%.3f")
    if cascade:
        # identical initial state: burned reservoir + Idukki slab pool
        # (LISFLOOD cannot drain the pool through the breach - noted)
        write_asc(out / "water.asc", dom.depth0, tr)

    # Mullaperiyar hydrograph -> per-point m^2/s series
    t, q = np.load(ROOT / "outputs" / "baseline_142" / "hydrograph.npy")
    t5 = np.arange(0, 86400 + 1, 300.0)
    q5 = np.interp(t5, t, q)
    pts = list(zip(dom.inj_rows.tolist(), dom.inj_cols.tolist()))
    npts = len(pts)

    with open(out / "mulla.bci", "w") as f:
        for r, c in pts:
            x = tr.c + (c + 0.5) * tr.a
            y = tr.f + (r + 0.5) * tr.e
            f.write(f"P {x:.1f} {y:.1f} QVAR inflow\n")
        if cascade:
            for r, c in zip(dom.cheru_rows.tolist(),
                            dom.cheru_cols.tolist()):
                x = tr.c + (c + 0.5) * tr.a
                y = tr.f + (r + 0.5) * tr.e
                f.write(f"P {x:.1f} {y:.1f} QVAR cheruthoni\n")

    with open(out / "mulla.bdy", "w") as f:
        f.write("Mullaperiyar Froehlich 142ft breach hydrograph\n")
        f.write("inflow\n")
        f.write(f"{len(t5)} seconds\n")
        for ti, qi in zip(t5, q5):
            f.write(f"{qi / (tr.a * npts):.6f} {ti:.0f}\n")
        if cascade:
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            import breach
            tc, qc, qcp, _ = breach.hydrograph(breach.cheruthoni_spec())
            nch = len(dom.cheru_rows)
            # zero until the prescribed trigger, then the full hydrograph
            tt = np.arange(0, 86400 + 1, 300.0)
            qq = np.interp(tt - trigger_t, tc, qc, left=0.0, right=0.0)
            f.write("cheruthoni\n")
            f.write(f"{len(tt)} seconds\n")
            for ti, qi in zip(tt, qq):
                f.write(f"{qi / (tr.a * nch):.6f} {ti:.0f}\n")

    # stage gauges
    g = []
    r, c = int(dom.inj_rows[0]), int(dom.inj_cols[0])
    g.append(("dam_toe", tr.c + (c + 0.5) * tr.a, tr.f + (r + 0.5) * tr.e))
    towns = ["vandiperiyar"] + (
        ["neriamangalam", "kalady", "aluva", "varappuzha"] if cascade else [])
    for name in towns:
        r, c = dom.gauges[name]
        g.append((name, tr.c + (c + 0.5) * tr.a, tr.f + (r + 0.5) * tr.e))
    rc = dom.ll_to_rc(9.820, 76.940)
    g.append(("idukki_pool", tr.c + (rc[1] + 0.5) * tr.a,
              tr.f + (rc[0] + 0.5) * tr.e))
    with open(out / "mulla.stage", "w") as f:
        f.write(f"{len(g)}\n")
        for _, x, y in g:
            f.write(f"{x:.1f} {y:.1f}\n")
    (out / "gauge_names.json").write_text(
        json.dumps([name for name, _, _ in g]))

    par = f"""# Mullaperiyar {args.scenario} - LISFLOOD-FP 8 (local inertial)
DEMfile        dem.asc
resroot        mulla
dirroot        results
sim_time       86400.0
initial_tstep  1.0
massint        60.0
saveint        3600.0
manningfile    manning.asc
bcifile        mulla.bci
bdyfile        mulla.bdy
stagefile      mulla.stage
acceleration
depththresh    0.01
"""
    if cascade:
        par += "startfile      water.asc\n"
    (out / "mulla.par").write_text(par)
    print(f"lisflood {args.scenario} project written to {out} "
          f"({npts} inflow points)")


if __name__ == "__main__":
    main()
