# Mullaperiyar dam-break flood simulation

A 2D shallow-water dam-break simulation of the Mullaperiyar Dam (Idukki
district, Kerala) over real Copernicus GLO-30 terrain, built for
**visualization and intuition only**. It is not an engineering study, it
says nothing about whether the dam will actually fail, and every number in
it is an order-of-magnitude estimate at best.

## Scenarios

| scenario | pool | Mullaperiyar breach | Idukki treatment |
|---|---|---|---|
| `baseline_142` | 142 ft (380 Mm³) | Froehlich (2008), ~175 m wide over ~2.5 h | pool **rises physically** behind the sealed dam line, spillway closed; rise + peak inflow reported |
| `cascade_142` | 142 ft (380 Mm³) | Froehlich (2008) | Cheruthoni dam (proxy for the Idukki cluster, 1,460 Mm³) breaches when the surge arrives; pool drains through the breach; flood routed to the sea |
| `sudden_152` | 152 ft (443 Mm³) | instantaneous collapse of 1/3 of the crest | cascade on |

Outputs per scenario in `outputs/<scenario>/`: `max_depth.tif` /
`arrival_time.tif` (EPSG:4326), `max_depth.png` / `arrival.png`,
`animation.mp4` (hillshade + depth), `animation_3d.mp4` (PyVista orbit +
valley fly-through), `report.md` (peak discharge, per-town arrival times
and peak depths), `run_meta.json` / `results.npz` (raw series).

## How to run

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python numpy numba rasterio matplotlib \
    imageio imageio-ffmpeg scipy requests tqdm pyvista scikit-image
.venv/bin/python src/fetch_dem.py            # downloads + conditions DEM
.venv/bin/python tests/smoke_test.py         # solver sanity (must pass)
.venv/bin/python src/run_scenario.py baseline_142 --res 90
.venv/bin/python src/viz2d.py outputs/baseline_142
.venv/bin/python src/viz3d.py outputs/baseline_142   # optional 3D
```

24 h of simulated flood takes ~1–3 min per scenario on a 16-core laptop
(numba, ~50k adaptive steps on a 1408×1053 grid).

## Model summary

- 2D shallow-water equations, first-order finite volume, Rusanov flux with
  Audusse hydrostatic reconstruction (exactly well-balanced: a lake at rest
  stays at rest to machine precision — see `tests/smoke_test.py`).
- CFL 0.4 adaptive timestep, 1 cm wet/dry thin film, semi-implicit Manning
  friction (n = 0.06 above 300 m ASL — forested gorge — else 0.035).
- Breach outflow from a broad-crested weir over a time-growing breach
  (Froehlich 2008 geometry/formation time), coupled to a level-pool
  stage–storage curve calibrated to the scenario volume; the hydrograph is
  injected at the gorge floor below the dam. Computed peaks: **~50,500 m³/s**
  (Froehlich @142 ft), **~65,300 m³/s** (sudden 1/3-crest @152 ft),
  **~387,000 m³/s** (Cheruthoni cascade, K₀ = 1.3).
- Mass ledger closes to +0.000 % in every production run (threshold 2 %).

## Limitations — read before quoting any number

1. **The DEM is a 30 m surface model resampled to 90 m.** Copernicus GLO-30
   records canopy tops, not the ground: the Periyar's ~50–150 m wide
   forested canyon is smoothed into a 200–400 m wide, shallower valley.
   The simulated surge is therefore *wider, shallower and slower* than a
   real one. Published estimates put the surge at Idukki in well under an
   hour; this model needs ~5–7 h. Treat all arrival times as biased late —
   possibly by several ×.
2. **The DEM was hydrologically conditioned.** Depression filling with an
   ε-gradient (mean fill ≈ 2 m, locally tens of metres in the gorge) was
   required to stop canopy artifacts from ponding the entire flood. Real
   valley storage is consequently underestimated.
3. **The breach hydrograph is an empirical guess.** Froehlich (2008) is a
   regression over historical earthen-dam failures; Mullaperiyar is a
   composite gravity dam — failure width/timing could differ greatly. The
   "sudden" scenario brackets the fast end.
4. **The Idukki treatment is stylized.** The lake is given a nominal
   initial pool (surface 730 m ASL: the DSM-recorded 725 m water surface
   plus a 5 m slab standing in for unknown bathymetry — the true level on
   any given day is unknown), the three dams are modeled as solid 742 m
   wall blocks (verified watertight at build time), and the **spillway is
   assumed closed** — no gate operations, so the surge has nowhere to go
   and the pool simply rises. Volume must convey through the lake's
   narrow dendritic arms, so the mid-lake level lags the eastern arm by
   hours. The cascade trigger (1 Mm³ of arrival) and Cheruthoni's
   Froehlich failure are placeholders, not dam-safety analysis — the arch
   dam is designed for far larger loads. In cascade runs, downstream
   peaks are dominated by the Cheruthoni release, so `cascade_142` and
   `sudden_152` differ mainly in timing, not magnitude.
5. **No sediment, debris, dam-body obstruction, buildings, channel
   structures or spillway operations.** Bhoothathankettu barrage and every
   downstream regulator are transparent; bridges/embankments that would
   obstruct or redirect a real flood are absent from the DEM.
6. **Vembanad/coastal hydraulics are not modeled** (no tides, no proper
   estuary bathymetry — the DSM shows the water surface). Depths in the
   Aluva–Varappuzha–Kochi plain are especially soft; the backwater peak is
   still rising at the 24 h cutoff in the cascade scenarios (published
   cascade studies quote ~5 m at Varappuzha; this model reaches ~2–2.5 m
   and climbing at t = 24 h).
7. **Numerics:** first-order scheme on 90 m cells smears fronts; Manning n
   is a two-value map; the reservoir drains as a level pool; ~1–3 % of the
   released volume exits through a spurious spill path west of the Idukki
   sink mask (Kulamavu saddle in the smoothed DSM).

**None of this says anything about the actual probability or consequence
of a Mullaperiyar failure.** It is a physics-flavored animation of one
prescribed what-if, on smoothed terrain, with empirical breach parameters.

## Comparison with published studies

Two quantitative comparators: the 2011 IIT Roorkee dam-break analysis
(commissioned during the Kerala–Tamil Nadu dispute; figures via the
summary at expert-eyes.org) and George et al. (2022), *Dam Break Hazard
Mapping: A Case Study of Mullaperiyar Dam*, IRJET 9(7) — a 1D HEC-RAS
study on an ALOS DEM covering both the Mullaperiyar breach and an Idukki
cascade. A smaller-breach HEC-RAS study (IJRASET) is also cited where
relevant.

### Assumptions side by side

| assumption | IIT Roorkee (2011) | George et al., IRJET (2022) | IJRASET study | **this model (2026)** |
|---|---|---|---|---|
| Hydraulic model | HEC-RAS 4.1, 1D St. Venant | HEC-RAS, 1D St. Venant | HEC-RAS, 1D | 2D shallow-water FV (Rusanov + hydrostatic reconstruction), numba |
| Terrain | surveyed cross-sections | ALOS 30 m DEM | SRTM/HEC-GeoRAS | Copernicus GLO-30 canopy DSM, ε-pit-filled, 90 m cells |
| Pool at failure | FRL 152 ft | 152 ft; 2018-flood peak inflow added | 152 ft | 142 ft (380 Mm³) and 152 ft (443 Mm³) variants |
| Mullaperiyar breach | full development in **12 min** | overtopping, fast | conservative small breach | Froehlich 2008 (175 m, **151 min**) and instant 1/3-crest (122 m) |
| Idukki pool response | overtopping wave assessed | inflow summed into Idukki analysis | — | pool rises physically behind sealed dam line, **spillway closed** |
| Idukki/Cheruthoni breach | — | partial arch-dam breach | — | Froehlich on full live storage (1,460 Mm³), K₀ = 1.3 |
| Downstream hydraulics | to Idukki | to the Arabian Sea, 1D channel | 36 km reach | to the sea, 2D overland spread, no embankments/structures |
| Volume conservation | n/r | n/r | n/r | mass ledger closes to +0.000 % |

### Results side by side

| result | published | **this model** | agreement |
|---|---|---|---|
| Peak breach outflow, Mullaperiyar | 89,121 m³/s (IIT-R / IRJET); 15,405 m³/s (IJRASET) | 50,540 (Froehlich) / 65,287 (sudden) m³/s | inside the published range; ~30 % below the high end, which assumes a 12-min breach vs Froehlich's 151 min |
| Depth just below the dam | 45.3 m (IRJET); 40.3 m (IJRASET) | 39.7–40.6 m | **matches** |
| Arrival at Vandiperiyar (7.4 km) | 25 min (IRJET); Vallakkadavu (3.6 km) 26 min (IIT-R) | 35 min (sudden) / 111 min (baseline) | sudden ~1.4× later; baseline later mostly because the Froehlich breach opens slowly |
| Peak depth at Vandiperiyar | 28.4 m (IRJET) | 29.0–30.1 m | **matches (±5 %)** |
| Arrival at Idukki reservoir | 122 min (IRJET) / 128 min (IIT-R) | 311 min (sudden) / 410 min (baseline) | **2.4–3.4× slower** — the canopy-DSM canyon-widening bias of limitation 1, quantified |
| Share of released volume reaching Idukki | ~85 % (IRJET) | 78 % (297/380 Mm³ impounded) | matches |
| Idukki pool response | overtopping of Cheruthoni assessed (IIT-R) | +5.2 m rise (730 → ~735.2 m ASL, crosses FRL at ~15 h; spillway closed) | comparable premise: the surge exceeds FRL headroom |
| Idukki-breach peak (cascade) | 30,458 m³/s (IRJET, small partial breach) | 387,000 m³/s (Froehlich, full storage) | structurally different breach assumptions, opposite directions |
| Idukki breach → Aluva travel time | ~9.7 h (IRJET) | ~10.1 h (arrival 17.3 h − trigger 7.25 h) | **matches (±5 %)** |
| Lowland depths in the cascade | Aluva 16 m, Ernakulam 7.5 m (1D IRJET); ~5 m at Varappuzha (other cascade studies) | Neriamangalam 46.7 m; Kalady 3.2 m; Aluva 2.3 m; Varappuzha 1.9–2.4 m (still rising at 24 h) | gorge depths comparable; **plain depths 3–7× shallower** — 1D cross-sections confine the flow to the channel; this 2D model spreads it across the ~30 km Vembanad plain and sees no embankments |

The pattern: gorge-confined physics (depths, volume fractions, lowland
travel times) agrees well; quantities dominated by breach assumptions
(peak discharge, early arrivals) differ by roughly the ratio of assumed
breach-formation speeds; and the coastal-plain depths differ structurally
with 1D channel models — reality likely sits between the two. The
IIT Roorkee 12-minute breach at 30–45 km/h wave speed is itself an
aggressive worst-case assumption, not a measurement.

On Idukki in `baseline_142`, the pool rise is now *simulated* rather than
inferred: ~297 Mm³ impounds behind the closed dams and the mid-lake gauge
rises **~5.2 m** (nominal 730 → ~735 m ASL), crossing FRL 732.6 m around
hour 15 — in good agreement with the independent stage–storage arithmetic
(+6.3 m at full leveling for that volume). Started at FRL instead, the
same volume would sit ~6 m above FRL; with the modeled peak inflow of
~7,000–8,000 m³/s exceeding ordinary spillway practice for hours, holding
the level would be impossible — the premise the cascade scenarios make
explicit.

Sources: [IRJET 9(7) 2022 paper (PDF)](https://www.irjet.net/archives/V9/i7/IRJET-V9I792.pdf),
[IIT Roorkee summary](https://www.expert-eyes.org/mullaperiyar/dam_break_analysis.html),
[IJRASET HEC-RAS study](https://www.ijraset.com/research-paper/dam-break-analysis-of-mullaperiyar-dam-using-hec-ras),
[Idukki Reservoir Break Analysis](https://www.researchgate.net/publication/354293406_Idukki_Reservoir_Break_Analysis).

## HEC-RAS 2D cross-check (attempted, incomplete)

An attempt to run HEC-RAS 6.6 2D on the same terrain + hydrograph lives in
`src/export_hecras.py`, `src/hecras_quadmesh.py`, `src/hecras_inject_ec.py`,
`src/hecras_trim_orphans.py` (toolchain: `neeraip/hecras-v66-linux`,
gitignored under `vendor/`). Status:

- HEC's official Linux compute engines run natively here (bundled Muncie
  and BEFORE_RUN examples reproduce end-to-end).
- A fully headless project-authoring pipeline works: corridor mesh
  (structured quads validated against every RASMapper topology convention),
  30 m subgrid terrain, Manning land cover, breach hydrograph as a BC-line
  flow, boundary/event conditions injected into the plan HDF.
- One configuration simulated the full 24 h (99.7 % complete before an
  end-of-run thread crash), but with the inflow volume being deleted at
  the boundary; after fixing the face-orientation conventions the solver
  became unstable on the hand-authored geometry (run-to-run varying init
  segfaults - symptomatic of some remaining undocumented format
  expectation).
- Conclusion: headless authoring of RAS 2D geometry from scratch is ~90 %
  reproduced but not solver-stable. The reliable route is a one-time mesh
  export from RAS Mapper on Windows (Workflow A/B), which then re-runs
  natively on Linux - proven with the bundled examples.

## LISFLOOD-FP cross-check (completed)

Since headless HEC-RAS stalled, the independent-solver cross-check was
done with **LISFLOOD-FP 8.0.3** (Bristol/Sheffield, GPL, built from the
Zenodo source; local-inertial "acceleration" solver) on the **identical**
180 m conditioned terrain, identical Froehlich 142 ft hydrograph at the
same gorge cells, identical Manning map, dry Idukki start. Tooling:
`src/export_lisflood.py`, `src/lisflood_compare.py`; results in
`lisflood/comparison.md`.

| quantity | ours (full-SWE Rusanov FV) | LISFLOOD-FP (local inertial) |
|---|---|---|
| Vandiperiyar arrival | 1.87 h | 1.78 h |
| Vandiperiyar peak depth | 21.9 m | 20.1 m |
| Vandiperiyar time of peak | 4.7 h | 5.2 h |
| Idukki basin volume at 24 h | 208 Mm³ | 247 Mm³ |
| Volume conservation | +0.000 % | ~0.09 Mm³ cumulative drift |
| Wall clock (24 h sim) | 10 s | 93 s |

Two independently-developed engines with different governing-equation
approximations agree to **~5 % on arrival time and ~8 % on peak depth**
at the key gorge gauge, and to ~20 % on the volume delivered to Idukki
(the inertial solver routes slightly faster). This is strong evidence
that the numbers here are properties of the terrain + hydrograph, not
artifacts of our solver - while leaving limitation 1 (the canopy-DSM
terrain itself) fully in force for both.

## Repository layout

```
src/fetch_dem.py     DEM download, mosaic, UTM grids, ε-fill conditioning
src/breach.py        Froehlich/sudden breach hydrographs (level-pool weir ODE)
src/solver.py        2D SWE finite-volume solver (numba)
src/terrain.py       reservoir burn, dam plugs, sink mask, gauges, Manning map
src/run_scenario.py  scenario configs + CLI runner
src/viz2d.py         GeoTIFFs, map animation, static maps, report.md
src/viz3d.py         PyVista 3D orbit + fly-through animation
tests/smoke_test.py  well-balancedness / conservation / ledger tests
docs/PLAN.md         design decisions
```
