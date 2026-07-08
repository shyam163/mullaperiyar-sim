"""2D deliverables per scenario: GeoTIFFs, animation, report.

Color rules (dataviz): depth and arrival are magnitudes -> single-hue
sequential ramps (ColorBrewer Blues / Oranges, truncated so the faintest
visible value still reads against the grey hillshade). The depth norm is
fixed across every animation frame so color means the same thing at every
timestamp. Towns are direct-labeled in ink with a halo, never color-coded.
"""
import json
import sys
from pathlib import Path

import imageio.v2 as imageio
import matplotlib
import numpy as np
import rasterio
from rasterio.warp import Resampling, calculate_default_transform, reproject

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.colors import LightSource, ListedColormap, PowerNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import terrain  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ARRIVAL_THRESH = 0.1

# truncated sequential ramps: skip the near-white low end
BLUES = ListedColormap(plt.get_cmap("Blues")(np.linspace(0.35, 1.0, 192)))
ORANGES = ListedColormap(plt.get_cmap("Oranges")(np.linspace(0.30, 1.0, 192)))


def _dem_for(meta):
    return ROOT / "data" / f"dem_utm{meta['res']}.tif"


def _load(scen_dir: Path):
    meta = json.load(open(scen_dir / "scenario_meta.json"))
    res = np.load(scen_dir / "results.npz")
    with rasterio.open(_dem_for(meta)) as src:
        z = src.read(1)
        profile = src.profile.copy()
    return meta, res, z, profile


def write_geotiffs(scen_dir: Path):
    """max_depth.tif and arrival_time.tif in EPSG:4326."""
    meta, res, _, profile = _load(scen_dir)
    for name, data, resampling, nodata in [
        ("max_depth.tif", res["max_depth"].astype(np.float32),
         Resampling.bilinear, 0.0),
        ("arrival_time.tif",
         np.where(res["arrival"] >= 0, res["arrival"] / 3600.0,
                  -9999.0).astype(np.float32),
         Resampling.nearest, -9999.0),
    ]:
        src_crs = profile["crs"]
        src_tr = profile["transform"]
        dst_tr, w, h = calculate_default_transform(
            src_crs, "EPSG:4326", profile["width"], profile["height"],
            *rasterio.transform.array_bounds(profile["height"],
                                             profile["width"], src_tr))
        out = np.full((h, w), nodata, np.float32)
        reproject(data, out, src_transform=src_tr, src_crs=src_crs,
                  dst_transform=dst_tr, dst_crs="EPSG:4326",
                  resampling=resampling, src_nodata=nodata,
                  dst_nodata=nodata)
        prof = dict(driver="GTiff", height=h, width=w, count=1,
                    dtype="float32", crs="EPSG:4326", transform=dst_tr,
                    nodata=nodata, compress="deflate")
        with rasterio.open(scen_dir / name, "w", **prof) as dst:
            dst.write(out, 1)
    print(f"geotiffs written for {scen_dir.name}")


class MapFigure:
    """Hillshade base + labeled towns; depth layer updated per frame."""

    def __init__(self, z, dom, meta, vmax_depth):
        self.fig, self.ax = plt.subplots(
            figsize=(12.8, 9.6), dpi=100)
        ls = LightSource(azdeg=315, altdeg=45)
        hs = ls.hillshade(z, vert_exag=2.5, dx=dom.dx, dy=dom.dx)
        self.ax.imshow(hs, cmap="gray", vmin=0, vmax=1.15)
        self.norm = PowerNorm(gamma=0.5, vmin=ARRIVAL_THRESH,
                              vmax=vmax_depth)
        self.im = self.ax.imshow(
            np.full(z.shape, np.nan), cmap=BLUES, norm=self.norm,
            interpolation="nearest", zorder=2)
        cb = self.fig.colorbar(self.im, ax=self.ax, shrink=0.62, pad=0.01)
        cb.set_label("water depth [m]")
        # Idukki lake: light static tint marks the pre-existing water body
        # (the DSM surface is the pool datum); the simulated pool RISE
        # renders on top of it through the normal depth layer
        tint = np.zeros(z.shape + (4,), np.float32)
        tint[dom.sink_rows, dom.sink_cols] = (0.42, 0.60, 0.78, 0.55)
        self.ax.imshow(tint, interpolation="nearest", zorder=1)
        halo = [patheffects.withStroke(linewidth=2.5, foreground="black")]
        marks = dict(terrain.TOWNS)
        marks["mullaperiyar dam"] = terrain.DAM_MULLA
        marks["idukki/cheruthoni"] = terrain.DAM_CHERUTHONI
        for name, (lat, lon) in marks.items():
            r, c = dom.ll_to_rc(lat, lon)
            dam = "dam" in name
            self.ax.plot(c, r, marker="^" if dam else "o",
                         ms=7 if dam else 5, mfc="white", mec="black",
                         mew=0.8, zorder=5)
            self.ax.annotate(name, (c, r), xytext=(6, 5),
                             textcoords="offset points", fontsize=9,
                             color="white", path_effects=halo, zorder=6)
        self.stamp = self.ax.text(
            0.01, 0.99, "", transform=self.ax.transAxes, va="top",
            fontsize=14, color="white", path_effects=halo, family="monospace")
        self.title = self.ax.set_title(meta["scenario"], fontsize=13)
        self.ax.set_axis_off()
        self.fig.tight_layout(pad=0.5)

    def frame(self, depth, t_s):
        d = np.where(depth >= ARRIVAL_THRESH, depth, np.nan)
        self.im.set_data(d)
        hh, mm = int(t_s // 3600), int(t_s % 3600 // 60)
        self.stamp.set_text(f"t = {hh:02d}:{mm:02d}")
        self.fig.canvas.draw()
        buf = np.asarray(self.fig.canvas.buffer_rgba())[:, :, :3]
        return buf


def render_animation(scen_dir: Path, every=1, fps=18, vmax_depth=15.0):
    meta, res, z, _ = _load(scen_dir)
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(_dem_for(meta), meta["breach"]["volume"])
    snaps = sorted((scen_dir / "snapshots").glob("h_*.npz"))[::every]
    fig = MapFigure(z, dom, meta, vmax_depth)
    out = scen_dir / "animation.mp4"
    with imageio.get_writer(out, fps=fps, codec="libx264", quality=7,
                            macro_block_size=16) as w:
        for i, sp in enumerate(snaps):
            with np.load(sp) as f:
                depth = f["h_cm"].astype(np.float32) / 100.0
                t_s = float(f["t"])
            w.append_data(fig.frame(depth, t_s))
            if i % 50 == 0:
                print(f"  frame {i}/{len(snaps)}", flush=True)
    plt.close(fig.fig)
    print(f"wrote {out}")


def render_statics(scen_dir: Path, vmax_depth=15.0):
    """max_depth.png and arrival.png for the report."""
    meta, res, z, _ = _load(scen_dir)
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(_dem_for(meta), meta["breach"]["volume"])
    for key, cmap, vmax, label, fname in [
            ("max_depth", BLUES, vmax_depth, "max water depth [m]",
             "max_depth.png"),
            ("arrival", ORANGES, 24.0, "arrival time [h]", "arrival.png")]:
        fig = MapFigure(z, dom, meta, vmax_depth)
        data = res[key].astype(np.float32)
        if key == "arrival":
            data = np.where(data >= 0, data / 3600.0, np.nan)
            fig.im.set_cmap(cmap)
            fig.im.set_norm(matplotlib.colors.Normalize(0, vmax))
            fig.im.colorbar.set_label(label)
            fig.im.set_data(data)
        else:
            fig.im.set_data(np.where(data >= ARRIVAL_THRESH, data, np.nan))
        fig.stamp.set_text("")
        fig.fig.savefig(scen_dir / fname, dpi=100)
        plt.close(fig.fig)
    print(f"statics written for {scen_dir.name}")


def write_report(scen_dir: Path):
    meta, res, z, _ = _load(scen_dir)
    r = meta["result"]
    gauge_names = r["gauge_names"]
    gs = res["gauge_series"]
    sink_bins = res["sink_bins"]
    gdt = r["gauge_dt"]
    br = meta["breach"]

    lines = [
        f"# Scenario report: {meta['scenario']}",
        "",
        f"- Grid: {meta['res']} m | duration {r['duration']/3600:.0f} h | "
        f"{r['steps']} steps | wall {r['wall_s']:.0f} s",
        f"- Breach: {br['name']}, avg width {br['b_avg']:.0f} m, formation "
        f"{br['t_form']/60:.0f} min",
        f"- **Peak breach discharge {br['q_peak']:,.0f} m3/s** at "
        f"t = {br['t_peak']/60:.0f} min after failure start",
        f"- Reservoir volume released: {br['volume']/1e6:.0f} Mm3 "
        f"(pool at WSE {meta['wse']:.1f} m ASL)",
        f"- Mass ledger error: {r['ledger_error']*100:+.3f}% "
        + ("**[FLAG > 2%]**" if abs(r["ledger_error"]) > 0.02 else "(OK)"),
        "",
        "## Idukki reservoir",
    ]
    if r["sink_first_arrival_s"] is not None:
        pk_bin = int(np.argmax(sink_bins))
        lines += [
            "- Pool modeled as a physically rising basin behind a sealed "
            "rim; **spillway assumed closed** (no releases)",
            f"- Surge arrival at reservoir: "
            f"**t = {r['sink_first_arrival_s']/3600:.2f} h**",
            f"- Peak volume impounded: {r['sink_total']/1e6:.0f} Mm3",
            f"- Peak inflow: {sink_bins[pk_bin]/gdt:,.0f} m3/s "
            f"at t = {pk_bin*gdt/3600:.1f} h",
        ]
        if "idukki_pool" in gauge_names:
            kp = gauge_names.index("idukki_pool")
            dp = gs[:, kp]
            last = int(np.nonzero(dp > 0)[0].max()) if (dp > 0).any() else 0
            level0 = 725.0 + dp[0]   # initial pool surface (datum + slab)
            rise = dp - dp[0]        # rise above the initial pool
            if (rise > 0.05).any():
                pk = float(rise.max())
                lines.append(
                    f"- **Pool level rise: {pk:.1f} m** above the nominal "
                    f"{level0:.0f} m initial surface (peak level "
                    f"~{level0+pk:.1f} m ASL, FRL is 732.6 m) at "
                    f"t = {rise.argmax()*gdt/3600:.1f} h")
                marks = [f"{rise[min(hh*60, last)]:.1f}"
                         for hh in (6, 9, 12, 18, 24)]
                lines.append(
                    "- Rise at t = 6/9/12/18/24 h: "
                    + " / ".join(marks) + " m")
    else:
        lines.append("- Surge never reached the reservoir")
    casc = meta.get("cascade") or {}
    if casc:
        lines += [
            "",
            "## Cascade (Cheruthoni breach)",
            f"- Triggered at t = {casc['trigger_t']/3600:.2f} h",
            f"- Cheruthoni peak discharge {casc['q_peak']:,.0f} m3/s, "
            f"{casc['t_peak_after_trigger']/60:.0f} min after trigger",
        ]
    lines += ["", "## Towns", "",
              "| gauge | arrival (h) | peak depth (m) | time of peak (h) |",
              "|---|---|---|---|"]
    for k, name in enumerate(gauge_names):
        if name == "idukki_pool":
            continue   # reported in the Idukki section, not a town
        d = gs[:, k]
        if (d > ARRIVAL_THRESH).any():
            ta = np.argmax(d > ARRIVAL_THRESH) * gdt / 3600
            lines.append(f"| {name} | {ta:.2f} | {d.max():.2f} | "
                         f"{d.argmax()*gdt/3600:.1f} |")
        else:
            lines.append(f"| {name} | - | dry | - |")
    lines += [
        "",
        "_Arrival threshold 0.1 m. Gauges snapped to the lowest cell within "
        "300 m of the town coordinate; depths are at that cell, not a town "
        "average. Order-of-magnitude estimates only - see top-level "
        "README for limitations._",
    ]
    (scen_dir / "report.md").write_text("\n".join(lines))
    print(f"report written for {scen_dir.name}")


def main(scen_dir):
    scen_dir = Path(scen_dir)
    vmax = 25.0 if "cascade" in scen_dir.name or "sudden" in scen_dir.name \
        else 15.0
    write_geotiffs(scen_dir)
    render_statics(scen_dir, vmax_depth=vmax)
    write_report(scen_dir)
    render_animation(scen_dir, vmax_depth=vmax)


if __name__ == "__main__":
    main(sys.argv[1])
