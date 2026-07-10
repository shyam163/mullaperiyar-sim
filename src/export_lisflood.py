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
    OUT.mkdir(exist_ok=True)
    (OUT / "results").mkdir(exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(ROOT / "data" / f"dem_utm{RES}.tif", 380e6)
    tr = dom.transform
    assert terrain.LAKE_SLAB == 0.0, "run with MULLA_LAKE_SLAB=0"

    write_asc(OUT / "dem.asc", dom.z, tr)
    write_asc(OUT / "manning.asc", dom.n_map, tr, fmt="%.3f")

    # hydrograph -> per-point m^2/s series at the injection cells
    t, q = np.load(ROOT / "outputs" / "baseline_142" / "hydrograph.npy")
    t5 = np.arange(0, 86400 + 1, 300.0)
    q5 = np.interp(t5, t, q)
    pts = list(zip(dom.inj_rows.tolist(), dom.inj_cols.tolist()))
    npts = len(pts)

    with open(OUT / "mulla.bci", "w") as f:
        for k, (r, c) in enumerate(pts):
            x = tr.c + (c + 0.5) * tr.a
            y = tr.f + (r + 0.5) * tr.e
            f.write(f"P {x:.1f} {y:.1f} QVAR inflow\n")
    with open(OUT / "mulla.bdy", "w") as f:
        f.write("Mullaperiyar Froehlich 142ft breach hydrograph\n")
        f.write("inflow\n")
        f.write(f"{len(t5)} seconds\n")
        for ti, qi in zip(t5, q5):
            f.write(f"{qi / (tr.a * npts):.6f} {ti:.0f}\n")

    # stage gauges: dam_toe, vandiperiyar, mid-lake
    g = []
    r, c = int(dom.inj_rows[0]), int(dom.inj_cols[0])
    g.append(("dam_toe", tr.c + (c + 0.5) * tr.a, tr.f + (r + 0.5) * tr.e))
    r, c = dom.gauges["vandiperiyar"]
    g.append(("vandiperiyar", tr.c + (c + 0.5) * tr.a,
              tr.f + (r + 0.5) * tr.e))
    rc = dom.ll_to_rc(9.820, 76.940)
    g.append(("idukki_pool", tr.c + (rc[1] + 0.5) * tr.a,
              tr.f + (rc[0] + 0.5) * tr.e))
    with open(OUT / "mulla.stage", "w") as f:
        f.write(f"{len(g)}\n")
        for _, x, y in g:
            f.write(f"{x:.1f} {y:.1f}\n")
    (OUT / "gauge_names.json").write_text(
        json.dumps([name for name, _, _ in g]))

    par = f"""# Mullaperiyar baseline_142 - LISFLOOD-FP 8 (local inertial)
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
    (OUT / "mulla.par").write_text(par)
    print(f"lisflood project written to {OUT} ({npts} inflow points, "
          f"{len(t5)} hydrograph steps)")


if __name__ == "__main__":
    main()
