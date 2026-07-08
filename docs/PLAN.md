# Mullaperiyar Dam-Break Flood Simulation — Design & Plan

Status: user-provided spec is the approved design. This doc = distilled decisions + build order.
Purpose: visualization & intuition, order-of-magnitude only. NOT engineering-grade.

## Build order (user directive 2026-07-06)
solver + 2D outputs first and validated; 3D PyVista renders LAST; if PyVista has
headless/EGL trouble, ship 2D and leave `viz3d.py` ready to run — do not block.

## Key constants
- Dam: 9.5286N, 77.1394E. Height 53.6 m from foundation, crest length 365.7 m.
- Pool volumes: 380 Mm3 @142 ft (head 43.3 m), 443 Mm3 @152 ft (head 46.3 m). Base ~850 m ASL (verify from DEM).
- Domain bbox: lon 76.15–77.30, lat 9.40–10.25 → UTM 43N (EPSG:32643), 90 m production / 180 m coarse.
- Reservoir clip box: 9.52–9.58 N, 77.10–77.25 E; burn via flood-fill from seed near dam,
  WSE solved so DEM-integrated volume ≈ scenario volume (DSM shows flat water surface, so
  absolute "base+43" is unreliable — anchor on volume, not elevation).
- Idukki sink mask: flood fill z < ~735 m from seed ~(9.85N, 77.00E) in box 9.70–9.95N, 76.85–77.13E
  (Idukki FRL = 732.6 m; DSM shows reservoir surface ~730).
- Cheruthoni dam: 9.845N, 76.977E, H=138 m, live storage 1460 Mm3. Downstream = NORTH
  (verified on DSM: lake surface 725-733 m lies south of the dam line; the gorge at 550-690 m runs north).
- Towns (snap gauge to lowest cell within 300 m): Vandiperiyar 9.566/77.088, Neriamangalam 10.05/76.78,
  Kalady 10.17/76.44, Aluva 10.11/76.35, Varappuzha 10.07/76.27. Idukki entry = recorded dynamically
  (first sink-mask cell the surge reaches).

## Architecture (src/)
- `fetch_dem.py`  — tiles (cached data/tiles/) → mosaic → clip → EPSG:32643 @90 m & @180 m GeoTIFFs in data/.
- `breach.py`     — Froehlich 2008 (B=0.27·K0·V^0.32·H^0.04, tf=63.2·sqrt(V/(g·H²))), trapezoidal breach
                    growing linearly over tf; broad-crested weir Q=1.7·b·H^1.5 + side terms; V-shape
                    stage-storage fit (V=½·a·η²) calibrated to pool volume; forward-Euler ODE @1 s.
                    Sudden variant: rectangular 1/3-crest (122 m) full-depth at t=0. Same machinery for Cheruthoni.
- `solver.py`     — 2D SWE, first-order FV, Rusanov flux + Audusse hydrostatic reconstruction (well-balanced),
                    CFL 0.4, thin film 1 cm, semi-implicit Manning, two-pass flux/update numba prange kernels
                    (no write races), dry-dry face early-out, active-window bounding box. Mass ledger:
                    initial + injected − extracted(sink+reservoir drain) − boundary outflow − storage ≈ 0.
- `terrain.py`    — masks, dam wall raise, injection cells (lowest channel cells just downstream), Manning map
                    (n=0.06 where z>300 m else 0.035), gauges.
- `scenarios.py`  — the 3 scenario configs; run_scenario.py CLI.
- `viz2d.py`      — GeoTIFFs (reproject to EPSG:4326), hillshade+depth mp4, report.md.
- `viz3d.py`      — PyVista orbit + fly-through, 1080p, every 3rd snapshot. LAST.

## Physics decisions
- Hydrograph is authoritative for flood volume (injected at cells just downstream of dam);
  the burned reservoir is drained at the same rate for visual/mass consistency (matched volume).
  Dam line stays a wall in the DEM.
- Reservoir routing upstream of dam is NOT modeled (level-pool assumption in breach ODE).
- Idukki: baseline → sink (tally volume + peak inflow). Cascade → sink stays on for arriving surge
  (crude, disclosed) + Cheruthoni Froehlich hydrograph (overtopping K0=1.3) triggered when
  cumulative surge inflow > 1 Mm3, injected NORTH of Cheruthoni dam (downstream gorge).
- Boundaries: west edge transmissive (sea); others reflective (flood never reaches them).
- Snapshots: depth uint16 (cm) compressed npz every 300 s → 288 per run.

## Sanity anchors
- Peak breach Q in 50k–300k m3/s band (print peak + t_peak).
- Surge at Idukki in tens of minutes.
- Cascade scenarios ~5 m depth at Varappuzha.
- Mass balance error < 2% (flag if worse).

## Verification gates
1. Smoke test (synthetic valley): lake-at-rest exact, mass conserved, wet front OK.
2. Ultracode adversarial review of solver/breach/terrain before production runs.
3. Coarse 180 m run: flood must follow Periyar valley — no watershed jumps.
4. Production runs: ledger + anchors above.
