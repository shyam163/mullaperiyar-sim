# Mullaperiyar dam-break flood simulation

A 2D shallow-water dam-break simulation of the Mullaperiyar Dam (Idukki
district, Kerala) over real Copernicus GLO-30 terrain, built for
**visualization and intuition only**. It is not an engineering study, it
says nothing about whether the dam will actually fail, and every number in
it is an order-of-magnitude estimate at best.

## Scenarios

| scenario | pool | Mullaperiyar breach | Idukki treatment |
|---|---|---|---|
| `baseline_142` | 142 ft (380 Mm³) | Froehlich (2008), ~175 m wide over ~2.5 h | absorbs the surge (infinite sink; volume + peak inflow reported) |
| `cascade_142` | 142 ft (380 Mm³) | Froehlich (2008) | Cheruthoni dam (proxy for the Idukki cluster, 1,460 Mm³) breaches when the surge arrives; flood routed to the sea |
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
4. **The Idukki treatment is crude.** The reservoir is an infinite sink;
   the cascade trigger (1 Mm³ of surge arrival) and the assumption that
   Cheruthoni fails instantly-with-Froehlich-timing are placeholders, not
   dam-safety analysis. The Idukki pool level, spillway operations, and the
   arch dam's actual robustness (it is designed for far larger loads) are
   all ignored. In cascade runs, downstream peaks are dominated by the
   Cheruthoni release, so `cascade_142` and `sudden_152` differ mainly in
   timing, not magnitude.
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
