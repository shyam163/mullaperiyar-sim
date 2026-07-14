# Mullaperiyar dam-break simulation

2D shallow-water dam-break sim (Mullaperiyar Dam, Kerala) over Copernicus
GLO-30 terrain. **Visualization/intuition only — order of magnitude, not
engineering.** Project is COMPLETE and PUBLISHED; most future work will be
small edits to the blog/video assets or re-runs with different knobs.

**Read `docs/JOURNEY.md` before changing anything substantive** — it is
the full record: every design decision, bug, dead end, and verified number
(§1–15 simulation, §16 publication). `README.md` = public summary with
comparison tables. `docs/PLAN.md` = original design. `docs/video_script.md`
= the ~5:40 YouTube script.

## Environment & commands

- venv: `.venv/` = **uv-managed Python 3.12** (system 3.14 lacks numba).
  Always use `.venv/bin/python`.
- Reproduce: `src/fetch_dem.py` (DEM + ε-fill) → `tests/smoke_test.py`
  (must pass) → `src/run_scenario.py {baseline_142|cascade_142|sudden_152}
  --res 90` (~1–4 min each) → `src/viz2d.py outputs/<scen>` →
  `src/viz3d.py outputs/<scen>` (add `--region lat0,lat1,lon0,lon1
  --suffix _x` for full-res crops; `--flysec N` for a fly-through pass).
- Env knobs: `MULLA_LAKE_SLAB` (Idukki initial pool, default 5 m),
  `MULLA_BREACH_TF_MIN`, `MULLA_N_GORGE`, `MULLA_CARVE` +
  `MULLA_CARVE_MODE=grade|sill`.

## Non-obvious traps (cost hours the first time)

- Copernicus is a **canopy DSM**: the gorge is a staircase of fake pits;
  `fetch_dem.py` applies priority-flood ε-fill. Arrival times are still
  biased ~2.4× late vs surveyed-terrain studies — documented, not a bug.
- Idukki geography: lake is SOUTH of the dam line, downstream gorge NORTH.
  Dams are sealed with solid rectangular blocks (discs/rings leak);
  `terrain.build_domain` verifies watertightness with a 741.5 m fill.
- Idukki initial pool must be a FLAT surface (730 m), not a uniform-depth
  blanket (it avalanches). Spillway assumed closed per user.
- PyVista StructuredGrid point arrays are **Fortran order** — C-order
  scalars = corduroy striping.
- Mass ledger must stay +0.000 %; smoke tests are machine-precision.

## Publication surfaces (keep in sync when editing the post)

Three copies of the article: `blog_v2.html` (repo), the deploy bundle at
`/home/shyam/.claude/jobs/ba82049d/tmp/blogdeploy/site/mullaperiyar/index.html`
(recreate from blog_v2 + chrome if gone), and LIVE at
https://blog.quantumautomata.in/mullaperiyar/ (`blog.html` is the frozen v1).
Deploy = scp to ubuntu@152.67.163.191 (key
`/home/shyam/Documents/keys/ssh-key-2025-04-13.key`) then sudo mv into
`/var/www/blog.quantumautomata.in/` + chown www-data. ~14 other vhosts on
that box — never disturb. Cloudflare zone is SSL Full (strict): new
subdomains need `sudo certbot --nginx -d <sub>.quantumautomata.in`.
GitHub: `github.com/shyam163/mullaperiyar-sim` — push after committing.

## Verified headline numbers (do not re-derive casually)

baseline_142 (90 m): peak Q 50,540 m³/s; Vandiperiyar 1.85 h / 29 m;
Idukki arrival 6.8 h; 297 Mm³ impounded; pool +5.2 m (730→735.2, crosses
FRL 732.6 ~h15). Cascade: Cheruthoni 387,000 m³/s; Neriamangalam 9 h/47 m;
Varappuzha ~2 m still rising at 24 h. tf=12 min experiment: 77,397 m³/s,
Idukki 316 min. LISFLOOD-FP cross-check agrees within 5–8 % (gorge).
Worst-case Ernakulam corridor: median 1.8 / p95 3.1 / max 4.6 m;
Pathanamthitta: zero (different basin). Ledger +0.000 % everywhere.
