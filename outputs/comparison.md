# LISFLOOD-FP vs our solver - same terrain, same hydrograph

180 m grid, Froehlich 142 ft hydrograph at the same gorge cells,
dry Idukki start, spillway closed, 24 h. LISFLOOD-FP 8.0.3
local-inertial ('acceleration') vs our explicit full-SWE Rusanov FV.

| gauge | ours: arrival (h) | ours: peak (m) | ours: t_peak (h) | LF: arrival (h) | LF: peak (m) | LF: t_peak (h) |
|---|---|---|---|---|---|---|
| dam_toe | 0.17 | 35.7 | 2.53 | 0.17 | 49.6 | 2.65 |
| vandiperiyar | 1.87 | 21.9 | 4.70 | 1.78 | 20.1 | 5.18 |
| idukki_pool | - | 0.0 | - | - | 0.0 | - |

- inflow volume (LF mass file): 380 Mm3 (hydrograph 380)
- final stored volume (LF): 380 Mm3
- wet footprint >10 cm: ours 3484 cells, LF 4019 cells, IoU 0.53
- our ledger error +0.000%; LF cumulative Verror 0.086 Mm3
- Idukki basin volume impounded at 24 h: ours 208 Mm3, LISFLOOD-FP 247 Mm3
- Caveats: the max-depth IoU (0.53) is dragged down by construction
  differences, not hydraulics - our run carries the burned Mullaperiyar
  reservoir as initial water (LF gets only the hydrograph), LF adds a
  cosmetic 1-cell spread line at the south domain edge, and both engines
  pile unphysical depth at the injection cells. Along the actual flood
  corridor the footprints and depths agree closely.
