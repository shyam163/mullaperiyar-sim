"""Write a structured-quad 2D mesh skeleton g01.hdf for HEC-RAS.

The hecras-v66-linux Workflow C (Voronoi from seeds) produces degenerate
geometry that segfaults RasUnsteady. This sidesteps it: uniform square
cells over the corridor polygon (staircase boundary, like our own FV
solver), written as the mesh-topology datasets that ras_preprocess's
verified Workflow B expects in an existing g01.hdf. The preprocessor then
computes the hydraulic tables on top.

Conventions (matched to RASMapper HDFs as read by read_mesh_from_g01hdf):
  - interior cells first (same order as Cell Points), ghost cells after,
    one ghost per perimeter face;
  - ghost cells have exactly one face and two face points;
  - FacePoints Is Perimeter: -1 for perimeter FPs, 0 otherwise.

Usage: python src/hecras_quadmesh.py [--res 180]
Reads hecras/project/{mulla.g01 [perimeter], gauges.json} written by
export_hecras.py, writes hecras/project/mulla.g01.hdf.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT / "hecras" / "project"
NAME = "mulla"
AREA = "Perimeter 1"


def read_perimeter_and_seeds():
    """Parse perimeter + seeds back out of the .g01 text we generated."""
    lines = (PROJ / f"{NAME}.g01").read_text().splitlines()
    peri, seeds, i = [], [], 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("Storage Area Surface Line="):
            n = int(ln.split("=")[1].strip())
            for j in range(1, n + 1):
                raw = lines[i + j]
                peri.append((float(raw[0:16]), float(raw[16:32])))
            i += n
        elif ln.startswith("Storage Area 2D Points="):
            n = int(ln.split("=")[1].strip())
            j = i + 1
            got = 0
            while got < n:
                raw = lines[j]
                for k in range(0, len(raw) - 31, 32):
                    seeds.append((float(raw[k:k + 16]),
                                  float(raw[k + 16:k + 32])))
                    got += 1
                    if got == n:
                        break
                j += 1
            i = j
        i += 1
    return np.array(peri), np.array(seeds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, default=180)
    a = ap.parse_args()
    S = float(a.res)

    peri, seeds = read_perimeter_and_seeds()
    print(f"perimeter {len(peri)} pts, seeds {len(seeds)}")

    # snap seeds onto integer grid indices; dedupe so the mesh spacing S
    # may be coarser than the seed spacing in the .g01
    x0 = seeds[:, 0].min() - S / 2
    y0 = seeds[:, 1].min() - S / 2
    ci0 = np.floor((seeds[:, 0] - x0) / S).astype(np.int64)
    cj0 = np.floor((seeds[:, 1] - y0) / S).astype(np.int64)
    uniq = np.unique(np.c_[ci0, cj0], axis=0)
    ci, cj = uniq[:, 0], uniq[:, 1]
    seeds = np.c_[x0 + (ci + 0.5) * S, y0 + (cj + 0.5) * S]
    NI, NJ = ci.max() + 1, cj.max() + 1
    occ = np.full((NI, NJ), -1, np.int64)     # cell grid -> interior index
    occ[ci, cj] = np.arange(len(seeds))
    centers = seeds.copy()
    print(f"{len(uniq)} cells after snapping to {S:.0f} m grid")

    # ---- corners (face points): grid nodes touching any occupied cell
    corner_id = np.full((NI + 1, NJ + 1), -1, np.int64)
    occb = occ >= 0
    touch = np.zeros((NI + 1, NJ + 1), bool)
    touch[:-1, :-1] |= occb
    touch[1:, :-1] |= occb
    touch[:-1, 1:] |= occb
    touch[1:, 1:] |= occb
    idx = np.nonzero(touch)
    corner_id[idx] = np.arange(len(idx[0]))
    n_fp = len(idx[0])
    fp_xy = np.c_[x0 + idx[0] * S, y0 + idx[1] * S]

    # ---- faces: x-faces (normal +x) between (i-1,j) and (i,j) at x=i*S,
    #             y-faces (normal +y) between (i,j-1) and (i,j) at y=j*S
    faces = []      # (cellL, cellR, fpA, fpB, nx, ny) - cellR = +normal side
    ghosts = []     # ghost centers
    n_int = len(seeds)

    def ghost(center):
        ghosts.append(center)
        return n_int + len(ghosts) - 1

    # RASMapper conventions (verified on the bundled examples):
    #   1. face normal = clockwise rotation of (fp1 - fp0)
    #   2. normal points from cell0 toward cell1
    #   3. ghost cells are always cell1 (normal therefore points OUTWARD
    #      at the perimeter) - violate this and inflow BCs pour water
    #      into ghost-land, where the solver deletes it
    for axis in (0, 1):
        if axis == 0:   # x-faces at x=i*S between (i-1,j) and (i,j)
            A = np.pad(occ, ((1, 1), (0, 0)), constant_values=-1)
            for i in range(NI + 1):
                left = A[i, :]
                right = A[i + 1, :]
                js = np.nonzero((left >= 0) | (right >= 0))[0]
                for j in js:
                    L, R = int(left[j]), int(right[j])
                    lo, hi = corner_id[i, j], corner_id[i, j + 1]
                    if L >= 0 and R >= 0:
                        # tangent +y -> normal +x = L->R
                        faces.append((L, R, lo, hi, 1.0, 0.0))
                    elif L < 0:   # ghost on -x side: normal -x, fps top->bottom
                        gh = ghost(centers[R] - [S, 0])
                        faces.append((R, gh, hi, lo, -1.0, 0.0))
                    else:         # ghost on +x side: normal +x
                        gh = ghost(centers[L] + [S, 0])
                        faces.append((L, gh, lo, hi, 1.0, 0.0))
        else:          # y-faces at y=j*S between (i,j-1) and (i,j)
            A = np.pad(occ, ((0, 0), (1, 1)), constant_values=-1)
            for j in range(NJ + 1):
                lo_ = A[:, j]
                hi_ = A[:, j + 1]
                is_ = np.nonzero((lo_ >= 0) | (hi_ >= 0))[0]
                for i in is_:
                    L, R = int(lo_[i]), int(hi_[i])
                    ca, cb = corner_id[i, j], corner_id[i + 1, j]
                    if L >= 0 and R >= 0:
                        # normal +y needs tangent -x: fps right->left
                        faces.append((L, R, cb, ca, 0.0, 1.0))
                    elif L < 0:   # ghost on -y side: normal -y, tangent +x
                        gh = ghost(centers[R] - [0, S])
                        faces.append((R, gh, ca, cb, 0.0, -1.0))
                    else:         # ghost on +y side: normal +y
                        gh = ghost(centers[L] + [0, S])
                        faces.append((L, gh, cb, ca, 0.0, 1.0))

    faces = np.array(faces, dtype=object)
    n_faces = len(faces)
    n_cells = n_int + len(ghosts)
    print(f"quad mesh: {n_int} interior + {len(ghosts)} ghost cells, "
          f"{n_faces} faces, {n_fp} face points")

    fci = np.array([[f[0], f[1]] for f in faces], np.int32)
    ffp = np.array([[f[2], f[3]] for f in faces], np.int32)
    fnl = np.array([[f[4], f[5], S] for f in faces], np.float32)

    # ---- per-cell face-point rings (CCW) and face lists
    cfp = np.full((n_cells, 5), -1, np.int32)
    for k in range(n_int):
        i, j = ci[k], cj[k]
        ring = [corner_id[i, j], corner_id[i + 1, j],
                corner_id[i + 1, j + 1], corner_id[i, j + 1]]
        cfp[k, :4] = ring
    ghost_face = {}
    for fi, f in enumerate(faces):
        for cell in (f[0], f[1]):
            if cell >= n_int:
                ghost_face[cell] = fi
    for cell, fi in ghost_face.items():
        cfp[cell, :2] = ffp[fi]

    # perimeter face points: those of faces touching a ghost
    fp_perim = np.zeros(n_fp, np.int32)
    gmask = (fci >= n_int).any(axis=1)
    fp_perim[np.unique(ffp[gmask])] = -1

    # Cells Face and Orientation Info/Values (read but unused by the
    # preprocessor's reader; still write a consistent version)
    cell_faces = [[] for _ in range(n_cells)]
    for fi in range(n_faces):
        cell_faces[fci[fi, 0]].append((fi, 1))
        cell_faces[fci[fi, 1]].append((fi, 0))
    info = np.zeros((n_cells, 2), np.int32)
    vals = []
    off = 0
    for c in range(n_cells):
        info[c] = (off, len(cell_faces[c]))
        vals.extend(cell_faces[c])
        off += len(cell_faces[c])
    vals = np.array(vals, np.int32)

    all_centers = np.vstack([centers] + ([np.array(ghosts)] if ghosts
                                         else []))

    # ---- write skeleton g01.hdf
    import rasterio
    wkt = rasterio.crs.CRS.from_epsg(32643).to_wkt()
    out = PROJ / f"{NAME}.g01.hdf"
    with h5py.File(out, "w") as f:
        f.attrs["Projection"] = np.bytes_(wkt)
        f.attrs["Units System"] = np.bytes_("US Customary")
        geo = f.create_group("Geometry")
        geo.attrs["Title"] = np.bytes_(NAME)
        geo.attrs["SI Units"] = np.bytes_("False")
        fa = geo.create_group("2D Flow Areas")
        att_dt = [("Name", "S16"), ("Mann", "<f4"),
                  ("Spacing dx", "<f4"), ("Spacing dy", "<f4"),
                  ("Cell Count", "<i4")]
        fa.create_dataset("Attributes", data=np.array(
            [(AREA.encode(), 0.0892, S, S, n_int)], att_dt))
        fa.create_dataset("Cell Points", data=seeds.astype(np.float64))
        ag = fa.create_group(AREA)
        ag.create_dataset("Perimeter", data=peri.astype(np.float64))
        ag.create_dataset("Cells Center Coordinate",
                          data=all_centers.astype(np.float64))
        ag.create_dataset("Cells Face and Orientation Info", data=info)
        ag.create_dataset("Cells Face and Orientation Values", data=vals)
        ag.create_dataset("Cells FacePoint Indexes", data=cfp)
        ag.create_dataset("Faces Cell Indexes", data=fci)
        ag.create_dataset("Faces FacePoint Indexes", data=ffp)
        ag.create_dataset("Faces NormalUnitVector and Length", data=fnl)
        ag.create_dataset("FacePoints Coordinate",
                          data=fp_xy.astype(np.float64))
        ag.create_dataset("FacePoints Is Perimeter", data=fp_perim)

        # fp -> faces (orientation +1 when the fp is the face's first
        # point) and fp -> cells adjacency, per RASMapper conventions
        fp_faces = [[] for _ in range(n_fp)]
        for fi in range(n_faces):
            fp_faces[ffp[fi, 0]].append((fi, 1))
            fp_faces[ffp[fi, 1]].append((fi, -1))
        fp_cells = [[] for _ in range(n_fp)]
        for c in range(n_cells):
            for fp in cfp[c]:
                if fp >= 0:
                    fp_cells[fp].append(c)
        def pack(lists, width):
            info = np.zeros((n_fp, 2), np.int32)
            vals = []
            off = 0
            for i, lst in enumerate(lists):
                info[i] = (off, len(lst))
                vals.extend(lst)
                off += len(lst)
            return info, np.array(vals, np.int32).reshape(-1, width)
        i1, v1 = pack(fp_faces, 2)
        ag.create_dataset("FacePoints Face and Orientation Info", data=i1)
        ag.create_dataset("FacePoints Face and Orientation Values", data=v1)
        i2, v2 = pack(fp_cells, 1)
        ag.create_dataset("FacePoints Cell Info", data=i2)
        ag.create_dataset("FacePoints Cell Index Values",
                          data=v2.ravel())

        # BC lines (RASMapper always writes this group; Workflow B expects
        # it and computes the External Faces itself from the polylines)
        import json
        gauges = json.loads((PROJ / "gauges.json").read_text())
        bc = geo.create_group("Boundary Condition Lines")
        attrs_dt = [("Name", "S32"), ("SA-2D", "S16"), ("Type", "S8"),
                    ("Length", "<f4")]
        attrs, pinfo, pparts, ppts = [], [], [], []
        for lid, (name, key) in enumerate(
                [("Inflow", "_bc_inflow"), ("Outlet", "_bc_outlet")]):
            pa, pb = gauges[key]
            ln = float(np.hypot(pb[0] - pa[0], pb[1] - pa[1]))
            attrs.append((name.encode(), AREA.encode(), b"External", ln))
            pinfo.append([len(ppts), 2, lid, 1])
            pparts.append([len(ppts), 2])
            ppts.extend([pa, pb])
        bc.create_dataset("Attributes", data=np.array(attrs, attrs_dt))
        bc.create_dataset("Polyline Info", data=np.array(pinfo, np.int32))
        bc.create_dataset("Polyline Parts", data=np.array(pparts, np.int32))
        bc.create_dataset("Polyline Points",
                          data=np.array(ppts, np.float64))
    print("wrote", out)


if __name__ == "__main__":
    main()
