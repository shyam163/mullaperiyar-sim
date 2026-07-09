"""Remove orphan (zero-facepoint, zero-face) cells from a Workflow-C
HEC-RAS geometry HDF.

The hecras-v66-linux Voronoi builder leaves behind cells whose clipped
region vanished (interior seeds too close to the perimeter) and ghost
cells it created but never wired. Nothing references them, but
RasUnsteady's geometry reader segfaults on them. This filters every
cell-indexed dataset, remaps face->cell and facepoint->cell indexes, and
mirrors the fixed Geometry into the plan HDF.

Usage: python src/hecras_trim_orphans.py <workdir> [--name mulla]
"""
import argparse
from pathlib import Path

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--name", default="mulla")
    a = ap.parse_args()
    wd = Path(a.workdir)
    g01p = wd / f"{a.name}.g01.hdf"
    p01p = wd / f"{a.name}.p01.tmp.hdf"

    with h5py.File(g01p, "a") as f:
        fa = f["Geometry/2D Flow Areas"]
        area = fa["Attributes"][()]["Name"][0].decode().strip("\x00").strip()
        g = fa[area]
        cfp = g["Cells FacePoint Indexes"][()]
        n_cells = len(cfp)
        n_int = fa["Cell Points"].shape[0]
        keep = (cfp >= 0).any(axis=1)
        orph = np.nonzero(~keep)[0]
        if len(orph) == 0:
            print("no orphan cells; nothing to do")
            return
        keep_int = keep[:n_int]
        n_int_new = int(keep_int.sum())
        print(f"trimming {len(orph)} orphan cells "
              f"({int((~keep_int).sum())} interior) -> "
              f"{int(keep.sum())} cells, {n_int_new} interior")

        # old->new cell index map
        newidx = np.full(n_cells, -1, np.int64)
        newidx[keep] = np.arange(keep.sum())

        def filt(dsname, grp=g, mask=keep):
            data = grp[dsname][()][mask]
            del grp[dsname]
            grp.create_dataset(dsname, data=data)

        for name in ["Cells Center Coordinate", "Cells Center Manning's n",
                     "Cells Face and Orientation Info",
                     "Cells FacePoint Indexes", "Cells Minimum Elevation",
                     "Cells Surface Area", "Cells Volume Elevation Info"]:
            filt(name)
        for sub in ["Infiltration", "Percent Impervious"]:
            if sub in g:
                for name in list(g[sub]):
                    d = g[sub][name]
                    if isinstance(d, h5py.Dataset) and d.shape[:1] == (n_int,):
                        filt(name, grp=g[sub], mask=keep_int)

        # remap references
        fci = g["Faces Cell Indexes"][()]
        assert not np.isin(fci, orph).any(), "face references an orphan"
        data = newidx[fci].astype(np.int32)
        del g["Faces Cell Indexes"]
        g.create_dataset("Faces Cell Indexes", data=data)

        fpcv = g["FacePoints Cell Index Values"][()]
        valid = fpcv >= 0
        assert not np.isin(fpcv[valid], orph).any(), "fp references orphan"
        fpcv[valid] = newidx[fpcv[valid]]
        del g["FacePoints Cell Index Values"]
        g.create_dataset("FacePoints Cell Index Values",
                         data=fpcv.astype(np.int32))

        # interior seed list + counts
        filt("Cell Points", grp=fa, mask=keep_int)
        att = fa["Attributes"][()]
        if "Cell Count" in att.dtype.names:
            att["Cell Count"][0] = n_int_new
            del fa["Attributes"]
            fa.create_dataset("Attributes", data=att)
        if "Cell Info" in fa:
            ci = fa["Cell Info"][()]
            ci[0, 1] = n_int_new
            del fa["Cell Info"]
            fa.create_dataset("Cell Info", data=ci)

    # mirror the fixed geometry into the plan HDF
    with h5py.File(p01p, "a") as fp, h5py.File(g01p, "r") as fg:
        del fp["Geometry"]
        fg.copy("Geometry", fp)
    print("geometry trimmed and mirrored into plan HDF")


if __name__ == "__main__":
    main()
