"""Fetch Copernicus GLO-30 DEM tiles and build the simulation grids.

Pipeline: download 1x1-degree COG tiles from the AWS open-data bucket
(no auth), mosaic them, clip to the simulation bounding box, then
reproject to UTM zone 43N (EPSG:32643) at 90 m and 180 m cell size.

Outputs:
    data/dem_utm90.tif   production grid
    data/dem_utm180.tif  coarse validation grid
"""
from pathlib import Path

import numpy as np
import rasterio
import requests
from rasterio.merge import merge
from rasterio.warp import Resampling, calculate_default_transform, reproject
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parent.parent
TILE_DIR = ROOT / "data" / "tiles"
DATA_DIR = ROOT / "data"

# Simulation bounding box (WGS84): Mullaperiyar dam to the Arabian Sea.
BBOX = (76.15, 9.40, 77.30, 10.25)  # (lon_min, lat_min, lon_max, lat_max)
TILES = ["N09_00_E076_00", "N09_00_E077_00", "N10_00_E076_00", "N10_00_E077_00"]
S3 = "https://copernicus-dem-30m.s3.amazonaws.com"
DST_CRS = "EPSG:32643"  # UTM 43N


def download_tiles():
    TILE_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for t in TILES:
        name = f"Copernicus_DSM_COG_10_{t}_DEM"
        path = TILE_DIR / f"{name}.tif"
        if not path.exists() or path.stat().st_size == 0:
            url = f"{S3}/{name}/{name}.tif"
            print(f"downloading {url}")
            r = requests.get(url, timeout=300)
            r.raise_for_status()
            path.write_bytes(r.content)
        paths.append(path)
    return paths


def build_grid(mosaic_path: Path, res: float, out_path: Path):
    """Reproject the clipped WGS84 mosaic to UTM 43N at `res` metres."""
    with rasterio.open(mosaic_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, DST_CRS, src.width, src.height, *src.bounds, resolution=res
        )
        profile = src.profile.copy()
        profile.update(
            crs=DST_CRS, transform=transform, width=width, height=height,
            dtype="float32", compress="deflate", nodata=None, driver="GTiff",
        )
        dst_arr = np.empty((height, width), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst_arr,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=DST_CRS,
            resampling=Resampling.bilinear,
        )
        # Copernicus uses 0 for ocean; keep as-is (sea level datum).
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(dst_arr, 1)
    print(f"wrote {out_path}  {width}x{height} @ {res} m")


def main():
    paths = download_tiles()
    srcs = [rasterio.open(p) for p in paths]
    mosaic, transform = merge(srcs, bounds=BBOX)
    profile = srcs[0].profile.copy()
    for s in srcs:
        s.close()
    profile.update(
        transform=transform, width=mosaic.shape[2], height=mosaic.shape[1],
        dtype="float32", compress="deflate", driver="GTiff", nodata=None,
    )
    clip_path = DATA_DIR / "dem_wgs84_clip.tif"
    with rasterio.open(clip_path, "w", **profile) as dst:
        dst.write(mosaic[0].astype(np.float32), 1)
    print(f"wrote {clip_path}  {mosaic.shape[2]}x{mosaic.shape[1]}")

    build_grid(clip_path, 90.0, DATA_DIR / "dem_utm90.tif")
    build_grid(clip_path, 180.0, DATA_DIR / "dem_utm180.tif")


if __name__ == "__main__":
    main()
