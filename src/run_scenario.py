"""Run one scenario end to end: terrain + breach + solver.

Usage:
    python src/run_scenario.py baseline_142 [--res 90] [--hours 24] [--tag x]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import breach  # noqa: E402
import terrain  # noqa: E402
from solver import Simulation  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

SCENARIOS = {
    # pool level, Mullaperiyar breach mode, does Idukki/Cheruthoni cascade?
    "baseline_142": dict(pool_ft=142, mode="froehlich", cascade=False),
    "cascade_142": dict(pool_ft=142, mode="froehlich", cascade=True),
    "sudden_152": dict(pool_ft=152, mode="sudden", cascade=True),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scenario", choices=SCENARIOS)
    ap.add_argument("--res", type=int, default=90, choices=[90, 180])
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    cfg = SCENARIOS[args.scenario]

    out_name = args.scenario + (f"_{args.tag}" if args.tag else "")
    out_dir = ROOT / "outputs" / out_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- breach hydrograph (authoritative flood volume)
    spec = breach.mullaperiyar_spec(cfg["pool_ft"], cfg["mode"])
    t_h, q_h, q_peak, t_peak = breach.hydrograph(spec)
    print(f"breach {spec.name}: B_avg={spec.b_avg:.0f} m "
          f"t_form={spec.t_form/60:.0f} min Qpeak={q_peak:,.0f} m3/s "
          f"@ {t_peak/60:.0f} min")

    # ---- terrain (pool volume must match the hydrograph's)
    dem = ROOT / "data" / f"dem_utm{args.res}.tif"
    dom = terrain.build_domain(dem, spec.volume_m3)

    sim = Simulation(dom.z, dom.dx, dom.n_map, out_dir,
                     duration=args.hours * 3600.0)
    sim.set_initial_water(dom.depth0)
    sim.add_injection(dom.inj_rows, dom.inj_cols, t_h, q_h,
                      drain_cells=(dom.res_rows, dom.res_cols))
    # Idukki: physically pooling basin behind the sealed rim, spillway
    # assumed CLOSED (no releases); monitored for reporting + cascade
    sim.set_basin(dom.sink_rows, dom.sink_cols)
    sim.gauges = dict(dom.gauges)
    # extra gauge just below the dam for the report
    sim.gauges["dam_toe"] = (int(dom.inj_rows[0]), int(dom.inj_cols[0]))
    # pool gauge mid-lake: depth here IS the rise above the DSM lake
    # surface (725 m datum)
    lake_seed = terrain.Domain.ll_to_rc(dom, 9.820, 76.940)
    sim.gauges["idukki_pool"] = lake_seed

    cascade_info = {}
    if cfg["cascade"]:
        cspec = breach.cheruthoni_spec()
        tc, qc, qcp, tcp = breach.hydrograph(cspec)

        def trigger(t_arr):
            print(f"CASCADE: Cheruthoni breach triggered at "
                  f"t={t_arr/3600:.2f} h", flush=True)
            sim.add_injection(dom.cheru_rows, dom.cheru_cols,
                              tc + t_arr, qc,
                              drain_cells=(dom.sink_rows, dom.sink_cols))
            cascade_info["trigger_t"] = t_arr
            cascade_info["q_peak"] = qcp
            cascade_info["t_peak_after_trigger"] = tcp

        sim.on_basin_arrival = trigger

    meta = dict(
        scenario=args.scenario, res=args.res,
        breach=dict(name=spec.name, b_avg=spec.b_avg, t_form=spec.t_form,
                    q_peak=q_peak, t_peak=t_peak,
                    volume=spec.volume_m3, head=spec.head_m),
        wse=dom.wse,
        gauges={g: dom.rc_to_ll(*rc) for g, rc in sim.gauges.items()},
    )
    result = sim.run()
    meta["result"] = result
    meta["cascade"] = cascade_info
    np.save(out_dir / "hydrograph.npy", np.stack([t_h, q_h]))
    with open(out_dir / "scenario_meta.json", "w") as f:
        json.dump(meta, f, indent=2, default=float)
    print("wrote", out_dir)


if __name__ == "__main__":
    main()
