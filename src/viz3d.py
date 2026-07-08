"""3D animation per scenario with PyVista (supplementary output).

Terrain: StructuredGrid from the DEM, 2.5x vertical exaggeration, colored
by a muted terrain colormap pre-blended with a hillshade (so no VTK light
setup is needed and the look survives any headless backend).
Water: the same grid warped to (terrain + depth), semi-transparent blue
ramp by depth, hidden where depth < 10 cm via NaN scalars.
Cameras: one full-domain orbit, then a fly-through following the Periyar
valley from the dam toward Kochi; both advance simulation time.

Usage: python src/viz3d.py outputs/<scenario> [--every 3]
If VTK cannot render offscreen on this machine, the script exits with a
clear message - ship the 2D outputs and run this later.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import terrain  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
VEXAG = 2.5
MIN_DEPTH = 0.1
DOWNSAMPLE = 2          # render grid = DEM grid / 2 (speed)
FPS = 12

# fly-through waypoints (lat, lon): dam -> Kochi backwaters
WAYPOINTS = [
    (9.5286, 77.1394), (9.566, 77.088), (9.700, 77.050), (9.800, 77.010),
    (9.845, 76.977), (9.930, 76.900), (10.05, 76.78), (10.11, 76.60),
    (10.17, 76.44), (10.11, 76.35), (10.07, 76.27), (10.05, 76.18),
]


def build_scene(dem_path, meta):
    import contextlib
    import io

    import matplotlib.pyplot as plt
    import pyvista as pv
    import rasterio
    from matplotlib.colors import LightSource

    with rasterio.open(dem_path) as src:
        z = src.read(1).astype(np.float64)
        tr = src.transform
    s = DOWNSAMPLE
    z = z[::s, ::s]
    ny, nx = z.shape
    dx = tr.a * s
    xs = tr.c + (np.arange(nx) + 0.5) * dx
    ys = tr.f - (np.arange(ny) + 0.5) * dx  # row 0 = north
    X, Y = np.meshgrid(xs, ys)

    with contextlib.redirect_stdout(io.StringIO()):
        dom = terrain.build_domain(dem_path, meta["breach"]["volume"])

    # --- terrain mesh with baked hillshade colors
    # NOTE: pv.StructuredGrid points iterate the meshgrid arrays in
    # FORTRAN order - every per-vertex array must be ravel(order="F")
    # or the colors/heights land on the wrong vertices (corduroy bands).
    from matplotlib.colors import LinearSegmentedColormap
    ls = LightSource(azdeg=315, altdeg=45)
    shade = ls.hillshade(z, vert_exag=2.0, dx=dx, dy=dx)[..., None]
    tcmap = LinearSegmentedColormap.from_list("earth", [
        "#42603f", "#6f8f4f", "#a8a35e", "#96825b", "#b5a284", "#e2ddd2"])
    tnorm = np.clip(z / 1500.0, 0, 1) ** 0.7
    rgb = tcmap(tnorm)[..., :3]
    rgb = np.clip(rgb * (0.35 + 0.85 * shade), 0, 1)
    terrain_grid = pv.StructuredGrid(X, Y, z * VEXAG)
    # rgb has shape (ny, nx, 3); each channel raveled F over (ny, nx)
    terrain_grid["rgb"] = np.stack(
        [(rgb[..., k] * 255).astype(np.uint8).ravel(order="F")
         for k in range(3)], axis=1)

    water_grid = pv.StructuredGrid(X, Y, z * VEXAG)
    water_grid["depth"] = np.full(z.size, np.nan)

    return dict(z=z, dx=dx, s=s, dom=dom, terrain=terrain_grid,
                water=water_grid, X=X, Y=Y)


def ll_to_xy(dom, lat, lon):
    from rasterio.warp import transform as wt
    xs, ys = wt("EPSG:4326", dom.crs, [lon], [lat])
    return xs[0], ys[0]


def smooth_path(pts, n):
    """Catmull-Rom-ish resample of waypoints to n points via cubic spline."""
    from scipy.interpolate import CubicSpline
    pts = np.asarray(pts)
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    cs = CubicSpline(d / d[-1], pts, axis=0)
    return cs(np.linspace(0, 1, n))


def render(scen_dir: Path, every=1, flysec=0):
    try:
        import pyvista as pv
    except ImportError:
        sys.exit("pyvista not installed - 3D render skipped")

    scen_dir = Path(scen_dir)
    meta = json.load(open(scen_dir / "scenario_meta.json"))
    dem_path = ROOT / "data" / f"dem_utm{meta['res']}.tif"
    sc = build_scene(dem_path, meta)
    z, dom, s = sc["z"], sc["dom"], sc["s"]

    snaps = sorted((scen_dir / "snapshots").glob("h_*.npz"))[::every]
    n = len(snaps)
    vmax = 25.0 if meta["scenario"] != "baseline_142" else 15.0

    try:
        pl = pv.Plotter(off_screen=True, window_size=(1920, 1088))
    except Exception as e:  # noqa: BLE001
        sys.exit(f"VTK offscreen init failed ({e}) - run later per README")

    pl.add_mesh(sc["terrain"], scalars="rgb", rgb=True, show_edges=False,
                smooth_shading=True)
    # truncated Blues: even 30 cm of water must read as water, not haze
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    blues = ListedColormap(plt.get_cmap("Blues")(np.linspace(0.45, 1.0, 160)))
    water_actor = pl.add_mesh(
        sc["water"], scalars="depth", cmap=blues, clim=(MIN_DEPTH, vmax),
        opacity=0.88, nan_opacity=0.0, show_edges=False, smooth_shading=True,
        show_scalar_bar=True,
        scalar_bar_args=dict(title="depth (m)", color="white",
                             position_x=0.90, position_y=0.30,
                             width=0.035, height=0.40, n_labels=3,
                             fmt="%.0f", title_font_size=20,
                             label_font_size=16))
    pl.set_background("#1a1f2b")

    # town labels
    marks = dict(terrain.TOWNS)
    marks["Mullaperiyar dam"] = terrain.DAM_MULLA
    marks["Idukki dams"] = terrain.DAM_CHERUTHONI
    lp, ln = [], []
    for name, (lat, lon) in marks.items():
        x, y = ll_to_xy(dom, lat, lon)
        r = np.clip(int((sc["Y"][0, 0] - y) / sc["dx"]), 0, z.shape[0] - 1)
        c = np.clip(int((x - sc["X"][0, 0]) / sc["dx"]), 0, z.shape[1] - 1)
        lp.append([x, y, z[r, c] * VEXAG + 900.0])
        ln.append(name)
    pl.add_point_labels(np.array(lp), ln, font_size=16, text_color="white",
                        shape_color="black", shape_opacity=0.45,
                        always_visible=True, show_points=True,
                        point_color="white", point_size=5)
    stamp = pl.add_text("t = 00:00", position="upper_left", font_size=14,
                        color="white")

    # --- camera paths
    center = np.array(sc["terrain"].center)
    ext = np.array([sc["X"].max() - sc["X"].min(),
                    sc["Y"].max() - sc["Y"].min()]).max()
    orbit_r = 0.62 * ext
    orbit_h = center[2] + 0.36 * ext
    ORBIT_ARC = np.radians(45.0)   # confined sweep, not a full circle
    ORBIT_BASE = np.pi             # centered due WEST: view from the sea

    # optional fly-through pass (flysec=0 -> arc only, the default): the
    # camera is interpolated smoothly while snapshots advance underneath
    n_fly = FPS * flysec
    fly_xy = smooth_path(
        [ll_to_xy(dom, la, lo) for la, lo in WAYPOINTS],
        n_fly) if n_fly > 0 else None
    if n_fly > 0:
        # two elevation series along the path: the valley floor (for the
        # focal point) and the highest terrain within ~3 km (for camera
        # clearance, so ridges beside the gorge cannot block the view)
        rad = max(int(3000.0 / sc["dx"]), 1)
        fly_floor, fly_clear = [], []
        for x, y in fly_xy:
            r = np.clip(int((sc["Y"][0, 0] - y) / sc["dx"]), 0,
                        z.shape[0] - 1)
            c = np.clip(int((x - sc["X"][0, 0]) / sc["dx"]), 0,
                        z.shape[1] - 1)
            fly_floor.append(z[r, c] * VEXAG)
            win = z[max(r - rad, 0):r + rad + 1,
                    max(c - rad, 0):c + rad + 1]
            fly_clear.append(win.max() * VEXAG)
        fly_floor = np.array(fly_floor)
        fly_clear = np.array(fly_clear)
        # smooth both so the camera does not bounce on gorge walls
        k = max(9, n_fly // 30)
        box = np.ones(2 * k + 1) / (2 * k + 1)
        fly_floor = np.convolve(np.pad(fly_floor, k, mode="edge"), box,
                                "same")[k:-k]
        # clearance must never smooth BELOW a ridge: running max first
        fly_clear = np.maximum(fly_clear, np.convolve(
            np.pad(fly_clear, k, mode="edge"), box, "same")[k:-k])

    import imageio.v2 as imageio
    out = scen_dir / "animation_3d.mp4"
    writer = imageio.get_writer(out, fps=FPS, codec="libx264", quality=7,
                                macro_block_size=16)

    def draw_frame(depth_cm, t_s):
        depth = depth_cm.astype(np.float32)[::s, ::s] / 100.0
        wsurf = np.where(depth >= MIN_DEPTH, (z + depth) * VEXAG + 2.0,
                         np.nan)
        pts = sc["water"].points.copy()
        pts[:, 2] = wsurf.ravel(order="F")
        # keep dry vertices on the terrain so triangles at the wet edge
        # do not stretch to NaN
        dryfix = np.isnan(pts[:, 2])
        pts[dryfix, 2] = (z * VEXAG).ravel(order="F")[dryfix] - 5.0
        sc["water"].points = pts
        dd = depth.copy()
        dd[depth < MIN_DEPTH] = np.nan
        sc["water"]["depth"] = dd.ravel(order="F")
        hh, mm = int(t_s // 3600), int(t_s % 3600 // 60)
        stamp.set_text("upper_left", f"t = {hh:02d}:{mm:02d}")
        pl.render()
        writer.append_data(pl.screenshot(return_img=True))

    # pass 1: 45-degree arc from the sea side while time advances
    for i, sp in enumerate(snaps):
        with np.load(sp) as f:
            hcm, t_s = f["h_cm"], float(f["t"])
        ang = ORBIT_BASE - ORBIT_ARC / 2 + ORBIT_ARC * i / max(n - 1, 1)
        pl.camera_position = [
            (center[0] + orbit_r * np.cos(ang),
             center[1] + orbit_r * np.sin(ang), orbit_h),
            center, (0, 0, 1)]
        draw_frame(hcm, t_s)
        if i % 48 == 0:
            print(f"arc {i}/{n}", flush=True)

    # pass 2: slow fly-through (n_fly frames); the camera glides while the
    # snapshot index advances underneath it
    ahead = max(8, n_fly // 36)
    last_si, hcm, t_s = -1, None, 0.0
    for i in range(n_fly):
        si = int(round(i * (n - 1) / max(n_fly - 1, 1)))
        if si != last_si:
            with np.load(snaps[si]) as f:
                hcm, t_s = f["h_cm"], float(f["t"])
            last_si = si
        j = min(i + ahead, n_fly - 1)
        cam = (fly_xy[i, 0], fly_xy[i, 1], fly_clear[i] + 5000.0)
        foc = (fly_xy[j, 0], fly_xy[j, 1], fly_floor[j] + 100.0)
        pl.camera_position = [cam, foc, (0, 0, 1)]
        draw_frame(hcm, t_s)
        if i % 60 == 0:
            print(f"fly {i}/{n_fly}", flush=True)

    writer.close()
    pl.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("scen_dir")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--flysec", type=int, default=0,
                    help="optional fly-through duration in seconds of "
                         "video (0 = sea-side arc only)")
    a = ap.parse_args()
    render(a.scen_dir, a.every, a.flysec)
