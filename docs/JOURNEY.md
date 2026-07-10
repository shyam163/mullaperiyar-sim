# Mullaperiyar dam-break simulation — complete technical documentation

Everything done in this project, in order, with every dead end, bug, fix
and number. Companion to `README.md` (user-facing summary) and
`docs/PLAN.md` (original design). Dates: 2026-07-06 → 2026-07-10.

---

## 1. Goal and ground rules

Build a 2D dam-break flood simulation of the Mullaperiyar Dam (9.5286 N,
77.1394 E, Idukki district, Kerala) over real topography, in Python, for
**visualization and intuition only** — explicitly not engineering-grade.
Domain: lon 76.15–77.30, lat 9.40–10.25 (dam → Arabian Sea, ~126 × 94 km).
Three scenarios; runtimes < 30 min each on a 16-core laptop; mass balance
logged and flagged over 2 %; published sanity anchors checked
(peak outflow 50–300k m³/s; surge at Idukki "in tens of minutes";
cascade ≈ 5 m at Varappuzha).

Physical constants used throughout:

| item | value |
|---|---|
| Dam height / crest length | 53.6 m / 365.7 m |
| Pool at 142 ft | 380 Mm³ (head 43.3 m) |
| Pool at 152 ft (FRL) | 443 Mm³ (head 46.3 m) |
| Cheruthoni dam (cascade proxy) | H = 138 m, live storage 1,460 Mm³, crest 650 m |
| Idukki reservoir FRL | 732.6 m ASL (DSM lake surface ≈ 725 m) |
| Manning n | 0.06 where z > 300 m, else 0.035 |
| Arrival threshold | 0.1 m depth |

## 2. Environment

- Nobara Linux (Fedora 43), 16 cores, 46 GB RAM, system Python 3.14.3.
- numba does not support 3.14 → `uv venv --python 3.12` (`.venv/`).
- Packages: numpy 2.4.6, numba 0.66 (16 threads), rasterio 1.5, scipy,
  matplotlib, imageio(-ffmpeg), pyvista 0.48, scikit-image, h5py, shapely,
  GDAL 3.13 (girder manylinux wheel — no system GDAL needed).

## 3. Terrain pipeline (`src/fetch_dem.py`)

1. Four Copernicus GLO-30 tiles (N09/N10 × E076/E077) fetched over plain
   HTTPS from the AWS open-data bucket (no auth), cached in `data/tiles/`.
2. Mosaic → clip to bbox → reproject EPSG:32643 (UTM 43N), bilinear, at
   90 m (1408×1053), 180 m (704×527), later 360 m.
3. **Hydrological conditioning — the first major discovery.** Copernicus
   is a *canopy DSM*: the forested, 50–150 m wide Periyar gorge is pinched
   shut into a staircase of closed depressions. First coarse run: the
   flood ponded 50 m deep below the dam and took >24 h to reach Idukki
   (published: ~2 h). Fixes, in sequence:
   - Plain morphological fill (`skimage.reconstruction`): ponds became
     exact-level flats → the wave *crawled*; Idukki arrival ~19 h.
   - **Priority-flood with ε-gradient** (Barnes 2014, `priority_flood_eps`,
     ε = 0.05 m/cell, applied only *inside* depressions so ocean/lakes stay
     untouched): 21.8 % of cells raised, mean 5.1 m. Idukki arrival → 9.4 h
     at 180 m, 6.8 h at 90 m. This is the conditioning kept for production.

## 4. Breach hydrographs (`src/breach.py`)

Level-pool ODE (Euler, 1 s) with broad-crested weir over a time-growing
breach; V-shaped stage–storage `V = ½aη²` calibrated to the pool volume.
Weir: `Q = 1.7·b·H^1.5 + 1.35·m·H^2.5` (SI). Froehlich (2008):
`B = 0.27·K₀·V^0.32·H^0.04`, `t_f = 63.2·√(V/(g·H²))`; side slope 1:1.
Sudden variant: 1/3 of crest (122 m), full depth at t = 0.

| spec | B_avg | t_form | Q_peak | t_peak | released |
|---|---|---|---|---|---|
| Mullaperiyar 142 ft Froehlich | 175 m | 151 min | 50,540 m³/s | 151 min | 380/380 Mm³ |
| Mullaperiyar 152 ft sudden | 122 m | 0 | 65,287 m³/s | 0 | 443/443 Mm³ |
| Mullaperiyar 142 ft, t_f forced 12 min | 175 m | 12 min | 77,397 m³/s | 12 min | 380/380 Mm³ |
| Cheruthoni cascade (K₀=1.3) | 366 m | 93 min | 387,111 m³/s | 93 min | 1,460/1,460 Mm³ |

Env knob: `MULLA_BREACH_TF_MIN` overrides Froehlich formation time.

## 5. Solver (`src/solver.py`)

First-order finite volume, 2D shallow-water:
- **Rusanov flux + Audusse et al. (2004) hydrostatic reconstruction** —
  provably well-balanced (lake-at-rest is an exact steady state).
- Adaptive dt, CFL 0.4; thin film H_EPS = 1 cm (mass kept, momentum
  zeroed); semi-implicit Manning friction `u/(1+dt·g n²|u|/h^{4/3})`;
  velocity safety cap 35 m/s; float64 state.
- numba `prange`, two passes/step (face fluxes → cell update) writing
  disjoint arrays; dry–dry faces early-out (cost ∝ wet area).
- **Active window**: bounding box of wet cells + margin, retightened every
  25 steps. Original version leaked mass at the window rim (CFL lets the
  front move 1 cell/step; the window refreshed every 25) → 50 % mass loss
  in the smoke test. Fix: grow window 1 cell/step between scans; airtight.
- **Boundaries**: 2-cell strips on all four domain edges are kept empty
  and swallowed into a `boundary` ledger bucket (kernel never updates the
  outermost cells → zero rim flux, every exit accounted).
- **Mass ledger**: initial + injected − drained − sink − boundary
  + clipped == storage; error printed hourly; +0.000 % in every
  production run. (Ledger sign on `clipped` was originally inverted —
  found by the adversarial review, then hand-verified; clipping was ~0 in
  all runs so no results were invalidated.)
- Injection: `add_injection(rows, cols, t, q, drain_cells=None)` — adds
  hydrograph mass at cells; optionally draws the same volume down from an
  upstream pool (Mullaperiyar reservoir for the main breach, Idukki pool
  for the cascade breach) so lakes visibly empty; hydrograph remains
  authoritative if the pool runs dry.
- Basin monitor (Idukki): tracks net volume gain vs the initial pool;
  records arrival (gain ≥ 10⁴ m³), fires the cascade callback at
  ≥ 1 Mm³; per-minute inflow bins.
- Gauges sampled every step into 60 s max-bins; max-depth and
  arrival-time rasters maintained incrementally in-kernel; snapshots
  every 300 s as compressed uint16-cm npz (289/run).

### Verification
- `tests/smoke_test.py`: lake-at-rest max|hu| 2.6e-14, surface flat to
  3.6e-15 m; closed-basin dam-break mass error 0.0; injection ledger
  1.8e-15. All pass.
- Analytic channel test: slope 0.003, n = 0.06, q = 9.26 m²/s → measured
  v = 2.36 m/s, h ≈ 3.5 m vs Manning normal 2.3 m/s, 4.0 m.
- Adversarial multi-agent review (5 lenses + 2 skeptics/finding):
  confirmed one stale doc statement (Cheruthoni downstream direction) and
  the `clipped` ledger sign; numba-race and geo-indexing lenses clean.

## 6. Domain assembly (`src/terrain.py`) — the geographic sagas

- **Reservoir burn**: flood-fill from a Thekkady seed inside a clip box
  (initially 9.52–9.58 N — clipped the real lake, extended to 9.44 N),
  with the WSE solved by bisection so the DEM-integrated volume matches
  the scenario pool (DSM shows the lake *surface*, so absolute
  base+43 m is meaningless). Production: WSE ≈ 885 m, 25.7 km², 380 Mm³.
  A 400 m disc raised to 930 m plugs the dam so the pool can't leak.
- **Injection cells**: lowest cells in a 350–900 m annulus NW of the dam
  (bearing 315°±90°), keeping only cells within 25 m of the channel floor.
- **Idukki mask**: naive fill below 736 m leaked to 278 km² (the DSM
  smooths the dam crests below pool level). Probing the DSM found: lake
  surface 724–733 m strictly *south* of the dam line; downstream gorge
  550–690 m runs *north* (the plan doc had this backwards); leak paths at
  the smoothed crests, an eastern corridor at 671–724 m, and the Kulamavu
  saddle (~734 m). Final mask: multi-seed fill from verified lake-surface
  points, post-filtered to z∈[695,736] and lat ≤ 9.848 → 50–56 km²
  (real FRL area ≈ 60 km²).
- **Idukki as a physically rising pool** (replaced the original
  infinite-sink treatment on user request; spillway assumed closed):
  1. Disc plugs at the three dams leaked. 2. An 810 m dilation-ring seal
  leaked through a NE corridor east of the lon-77° cut, then the pool
  drained into the gorge pocket *between* the dam line and the ring.
  3. Final: **solid rectangular wall blocks** across the Idukki-arch+
  Cheruthoni site and the Kulamavu saddle (crest 742 m), plus a build-time
  **verification flood-fill at 741.5 m** that must stay inside the
  reservoir footprint or `build_domain` raises.
  4. Initial pool: a *flat* surface at 730 m (`max(0, 730−z)`) — a
  uniform 5 m blanket was tried first and avalanched off the sloping
  shores (dt collapsed to 0.02 s). A dry-start lake was also tried: the
  surge then crawls across the flat lakebed (mid-lake gauge dry at 24 h)
  because a dry flat has no gravity-wave leveling; the nominal slab fixes
  this. `MULLA_LAKE_SLAB` env knob (production 5 m; cross-checks 0).
- **Cheruthoni cascade injection**: gorge floor *north* of the dam,
  1.4–2.8 km annulus (beyond the wall blocks), z ≈ 542 m at 90 m grid.
- **Synthetic channel knobs** (`MULLA_CARVE` depth, `MULLA_CARVE_MODE`
  grade|sill): steepest-descent thalweg from the injection cell to the
  basin (guaranteed to drain by the ε-fill), then either a linear graded
  bed to 726 m or a running-minimum sill-removal bed. See §11.
- Gauges snapped to the lowest cell within 300 m of each town.
  `MULLA_N_GORGE` overrides gorge Manning.

## 7. Scenarios and production results (90 m, 24 h, `src/run_scenario.py`)

| | baseline_142 | cascade_142 | sudden_152 |
|---|---|---|---|
| Mullaperiyar breach | Froehlich | Froehlich | instant 1/3 crest |
| peak Q | 50,540 m³/s | 50,540 | 65,287 |
| Idukki arrival | 6.8 h | 6.8 h | 5.1 h |
| Idukki response | +297 Mm³, pool +5.2 m (730→735.2, crosses FRL ~15 h) | trigger 7.25 h, pool drains via breach | trigger 5.57 h |
| Cheruthoni peak Q | — | 387,000 m³/s | 387,000 |
| Vandiperiyar | 1.85 h / 29.0 m | same | 0.58 h / 30.1 m |
| Neriamangalam | dry | 9.0 h / 46.7 m | 7.3 h / 46.7 m |
| Kalady | dry | 13.9 h / 3.2 m | 12.2 h / 3.2 m |
| Aluva | dry | 17.3 h / 2.3 m | 15.7 h / 2.3 m |
| Varappuzha | dry | 22.0 h / 1.9 m (rising at cutoff) | 20.3 h / 2.4 m |
| ledger error | +0.000 % | +0.000 % | +0.000 % |
| wall time | 59 s | 190 s | 219 s |

Downstream cascade peaks are identical between cascade_142 and sudden_152
because both are dominated by the same Cheruthoni release.

## 8. Comparison with published studies

Comparators: IIT Roorkee 2011 (HEC-RAS 1D, surveyed sections, 12-min
breach; via expert-eyes.org summary), George et al. 2022 IRJET (HEC-RAS
1D, ALOS 30 m), a small-breach IJRASET study. Highlights (full table in
README):
- Peak Q: published 15,405–89,121 m³/s; ours 50,540–65,287 — inside range.
- Depth below dam ~40 m and Vandiperiyar peak ~29 m: match within 5 %.
- Share of volume reaching Idukki: 85 % published vs 78–81 % ours.
- Idukki arrival: 122–128 min published vs 311–410 min ours (2.4–3.4×,
  see §11 for the decomposition).
- Idukki-breach travel time to Aluva: ~9.7 h vs ~10.1 h — match.
- Cascade lowland depths: 1D studies 7.5–16 m (channel-confined);
  ours 1.9–3.2 m (2D spreading over the ~30 km Vembanad plain, no
  embankments); the independent "~5 m at Varappuzha" anchor sits between.
- Post-hoc/simulated Idukki response: +5.2 m simulated rise agrees with
  independent stage–storage arithmetic (+6.3 m for that volume).

## 9. HEC-RAS 2D on Linux — attempted, ~90 % reproduced, not solver-stable

Toolchain: HEC's official Linux compute engines (RasGeomPreprocess,
RasUnsteady — run natively on Fedora, zero missing libs) + the
`neeraip/hecras-v66-linux` headless preprocessor. Muncie and BEFORE_RUN
examples reproduce end-to-end natively (gate passed). For our project:
- Fixed in the toolchain: `np.trapz` (numpy 2), a `%4d` fixed-width
  overflow in the .x01 writer at >9999 cells, stub .g01 BC-line parsers,
  orphan cells from its Voronoi mesher.
- Discoveries: mesh cells must be much coarser than the terrain raster
  (1 px/cell → degenerate volume curves → segfault; fixed with a 30 m
  terrain export); RASMapper mesh conventions (normal = CW-rotated
  tangent; normal points cell0→cell1; ghosts always cell1; CCW rings) —
  our replacement structured-quad mesh satisfies all of them, verified
  11,039/11,039 faces; Event Conditions in the plan HDF are what
  RasUnsteady actually reads (the .b01 hydrograph block is vestigial);
  the preprocessor hardcodes US units — exposed when the volume log
  reported exactly 380e6 ft³ as 8,723 acre-ft; countered with a
  consistent relabeling strategy (meters-as-feet, n × 1.486).
- Outcome: one configuration simulated 24 h to 99.7 % (inflow was being
  deleted at a mis-oriented boundary); after the orientation fixes the
  solver became non-deterministically unstable on the hand-authored
  geometry (init segfaults moving between runs) — the signature of some
  remaining undocumented dataset. Stopped per plan; scripts committed
  (`export_hecras.py`, `hecras_quadmesh.py`, `hecras_inject_ec.py`,
  `hecras_trim_orphans.py`). Reliable completion path: one-time RAS
  Mapper mesh export on Windows → proven Workflow A/B rerun on Linux.

## 10. LISFLOOD-FP 8.0.3 cross-check — completed

Built from the Zenodo source (cmake; 3-line `numa.h` shim + symlink to
`libnuma.so.1`; NetCDF off). Found and patched an **upstream bug**:
`LoadBCVar` stores `TimeSeries*` into `BCptr->PS_TimeSeries` while still
`push_back`-ing to the owning `std::vector` — any `.bdy` with >1 series
dangles the earlier pointers (segfault); fixed with `reserve(256)`.
Point-source units decoded from `iterateq.cpp`: QVAR values are m²/s
(per metre width; `H += q·dx·dt/dA`).

Same 180 m terrain, hydrograph, Manning, gauges as our matching runs:

| quantity | ours (full SWE) | LISFLOOD-FP (local inertial) |
|---|---|---|
| baseline: Vandiperiyar arrival / peak | 1.87 h / 21.9 m | 1.78 h / 20.1 m |
| baseline: Idukki basin volume @24 h | 208 Mm³ | 247 Mm³ |
| baseline: volume conservation | +0.000 % | ~0.09 Mm³ drift (380/380 in) |
| baseline: wall clock | 10 s | 93 s |
| cascade: Neriamangalam | 11.8 h / 50.5 m | 12.3 h / 36.3 m |
| cascade: Kalady | 16.5 h / 5.2 m | 19.3 h / 5.6 m |
| cascade: Aluva | 20.9 h / 2.0 m | 23.8 h / 1.0 m |
| cascade: footprint IoU | — | 0.72 |

Cascade forcing: Cheruthoni hydrograph prescribed at the trigger time
from our 180 m run (10.10 h); identical initial pools via `startfile`;
LISFLOOD cannot drain the pool through the breach (noted). The largest
depth disagreement (50 vs 36 m at Neriamangalam) is in the steep gorge
where the local-inertial approximation is known to diverge from full SWE.

## 11. Sensitivity experiments — decomposing the arrival-time gap

| run (90 m) | Vandiperiyar | Idukki | peak Q |
|---|---|---|---|
| baseline (t_f = 151 min) | 111 min | 410 min | 50,540 |
| t_f forced to 12 min | 41 min | 316 min | 77,397 |
| 12 min + graded 20 m channel | 45 min | 387 min | " |
| 12 min + sill-removal channel | 40 min | 303 min | " |
| sudden_152 + gorge n = 0.04 | 32 min | 278 min | 65,287 |
| published | 25 min | 122–128 min | 89,121 |

Lessons: breach time is ~60 % of the Idukki gap (85 % at Vandiperiyar);
roughness ~10 %; channel surgery backfires or barely helps — the first
carve *undercut the receiving lake by 9 m* (trench ponded at the mouth),
the graded rebuild was still slower than no channel because a uniform
0.0029 gradient is flatter than the natural steep reaches between sills
(and stores ~50 Mm³), and pure sill removal gains only 4 %. The residual
~2.4× is the canopy-DSM valley itself; closing it needs real DTM/
bathymetry data, not thalweg edits.

## 12. Visualization

- **2D (`src/viz2d.py`)**: per scenario — `max_depth.tif` /
  `arrival_time.tif` (EPSG:4326), `max_depth.png` / `arrival.png`,
  `animation.mp4` (hillshade + truncated-Blues depth ramp, fixed norm
  across frames, halo-labeled towns, timestamp), `report.md` (breach
  stats, Idukki pool section incl. spillway-closed statement and rise
  timeline, town table). Idukki lake drawn as a static tint *under* the
  dynamic layer.
- **3D (`src/viz3d.py`, PyVista off-screen)**: terrain StructuredGrid,
  2.5× exaggeration, hillshade pre-baked into vertex colors (headless-
  safe); water = same grid warped to surface, truncated Blues,
  NaN-hidden below 10 cm. Key gotcha: PyVista structured grids iterate
  points in **Fortran order** — C-order scalars produce corduroy
  striping. Camera evolved by user feedback: full orbit → terrain-
  clearance fly-through (max terrain within 2–3 km of path + altitude) →
  final form: **23 s arc confined to 45°, viewed from the Arabian Sea**
  looking east (fly-through now opt-in via `--flysec`).

## 13. Environment-variable knobs

| knob | effect |
|---|---|
| `MULLA_LAKE_SLAB` | Idukki initial pool depth datum (default 5 → surface 730 m; 0 = dry) |
| `MULLA_BREACH_TF_MIN` | override Froehlich formation time (minutes) |
| `MULLA_N_GORGE` | gorge Manning n (default 0.06) |
| `MULLA_CARVE` | synthetic channel depth (m); 0 = off |
| `MULLA_CARVE_MODE` | `grade` (linear to 726 m) or `sill` (running-min) |

## 14. Repository layout & reproduction

```
src/fetch_dem.py       tiles → UTM grids → ε-fill conditioning
src/breach.py          Froehlich / sudden hydrographs
src/solver.py          numba SWE solver
src/terrain.py         masks, walls, pools, gauges, carve knobs
src/run_scenario.py    scenario configs + CLI
src/viz2d.py, viz3d.py deliverable rendering
src/export_lisflood.py + lisflood_compare.py     LISFLOOD-FP cross-check
src/export_hecras.py + hecras_*.py               HEC-RAS attempt
tests/smoke_test.py    well-balancedness / conservation / ledger
docs/PLAN.md           original design   docs/JOURNEY.md  this file
outputs/<scenario>/    all deliverables  vendor/, lisflood/, hecras/  (gitignored)
```

Reproduce: `uv venv --python 3.12 .venv`; install deps (README);
`python src/fetch_dem.py`; `python tests/smoke_test.py`;
`python src/run_scenario.py <scenario> --res 90`;
`python src/viz2d.py outputs/<scenario>`; optional `viz3d.py`.

## 15. Honest-limitations register

Canopy DSM (arrival times biased late ~2.4× beyond breach assumptions);
ε-fill deletes real valley storage; no structures/embankments/spillway
operations anywhere; Idukki pool nominal (730 m) with spillway closed by
assumption; breach parameters are regressions over other dams; 1-cell
town gauges; no tides/estuary bathymetry (Vembanad soft; cascade peak
still rising at the 24 h cutoff); first-order scheme smears fronts;
~1–3 % of released volume exits a spurious Kulamavu-side path in early
runs (eliminated by the solid blocks); results are order-of-magnitude,
and none of this says anything about the probability of failure.
