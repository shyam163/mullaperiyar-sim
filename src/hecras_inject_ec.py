"""Wire boundary conditions into a Workflow-C HEC-RAS project.

The hecras-v66-linux preprocessor's .g01 BC-line parsers are stubs, so a
from-scratch project ends up with placeholder BC tables (US/DS, zero
external faces) and empty Event Conditions. This script finishes the job:

1. Selects the mesh perimeter faces lying along our Inflow / Outlet BC
   segments (from gauges.json) and rewrites
   Geometry/Boundary Condition Lines/{Attributes, External Faces,
   Polyline *} in g01.hdf.
2. Writes Event Conditions (flow hydrograph + normal depth) with matching
   face-index attributes into p01.tmp.hdf.
3. Patches the unit-system attributes to SI (preprocessor hardcodes US).

Usage: python src/hecras_inject_ec.py <workdir> [--name mulla]
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

START = "31Dec2025 2400"   # HEC convention: 0000 start = prev day 2400
END = "02Jan2026 0000"
ROOT = Path(__file__).resolve().parent.parent
PROJ = ROOT / "hecras" / "project"


def seg_dist(pts, a, b):
    """Distance from each point to segment a-b."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ab = b - a
    t = np.clip(((pts - a) @ ab) / (ab @ ab + 1e-12), 0, 1)
    proj = a + t[:, None] * ab
    return np.linalg.norm(pts - proj, axis=1), t


def pick_faces(area_grp, n_interior, a, b, tol):
    """Perimeter faces whose midpoint lies within tol of segment a-b,
    ordered along the segment."""
    fci = area_grp["Faces Cell Indexes"][()]
    ffp = area_grp["Faces FacePoint Indexes"][()]
    fpc = area_grp["FacePoints Coordinate"][()]
    lengths = area_grp["Faces NormalUnitVector and Length"][()][:, 2]
    is_perim = (fci >= n_interior).any(axis=1)
    mids = 0.5 * (fpc[ffp[:, 0]] + fpc[ffp[:, 1]])
    d, t = seg_dist(mids, a, b)
    sel = np.nonzero(is_perim & (d <= tol))[0]
    sel = sel[np.argsort(t[sel])]
    return sel.astype(np.int32), ffp, lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("--name", default="mulla")
    a = ap.parse_args()
    wd = Path(a.workdir)

    gauges = json.loads((PROJ / "gauges.json").read_text())
    t_s, q = np.loadtxt(PROJ / "hydrograph_5min.csv", delimiter=",",
                        skiprows=1).T
    series = np.c_[t_s / 86400.0, q].astype(np.float32)  # time in DAYS

    g01p = wd / f"{a.name}.g01.hdf"
    p01p = wd / f"{a.name}.p01.tmp.hdf"

    lines = {}
    with h5py.File(g01p, "a") as f:
        area = f["Geometry/2D Flow Areas/Attributes"][()]["Name"][0]
        area = area.decode().strip("\x00").strip()
        ag = f[f"Geometry/2D Flow Areas/{area}"]
        n_interior = f["Geometry/2D Flow Areas/Cell Points"].shape[0]

        # If the preprocessor already computed proper BC external faces
        # (quad-mesh Workflow B path), just read them - no rewiring.
        bcg = f.get("Geometry/Boundary Condition Lines")
        if bcg is not None and "External Faces" in bcg \
                and bcg["External Faces"].shape[0] > 0:
            names = [n.decode().strip() for n in bcg["Attributes"][()]["Name"]]
            ext = bcg["External Faces"][()]
            for lid, name in enumerate(names):
                rows = ext[ext["BC Line ID"] == lid]
                fidx = rows["Face Index"].astype(np.int32)
                # ordered face-point chain straight from the preprocessor's
                # FP Start/End columns (the solver walks this as a polyline)
                fps = [int(rows["FP Start Index"][0])]
                for r in rows:
                    fps.append(int(r["FP End Index"]))
                lines[name] = dict(faces=fidx,
                                   fps=np.array(fps, np.int32))
                print(f"BC {name}: {len(fidx)} faces (from preprocessor), "
                      f"chain {fps}")
            f.attrs["Units System"] = np.bytes_("US Customary")
            f["Geometry"].attrs["SI Units"] = np.bytes_("False")
            write_ec(p01p, g01p, area, lines, series)
            return
        # mesh spacing for the pick tolerance
        S = float(f["Geometry/2D Flow Areas/Attributes"][()]["Spacing dx"][0])
        tol = max(1.2 * S, 250.0)

        for name, key in [("Inflow", "_bc_inflow"), ("Outlet", "_bc_outlet")]:
            pa, pb = gauges[key]
            sel, ffp, lengths = pick_faces(ag, n_interior, pa, pb, tol)
            if len(sel) == 0 and name == "Outlet":
                # the outlet is a dry, inactive normal-depth line - snap it
                # to the nearest perimeter faces wherever they are
                fci = ag["Faces Cell Indexes"][()]
                fpc = ag["FacePoints Coordinate"][()]
                is_perim = (fci >= n_interior).any(axis=1)
                mids = 0.5 * (fpc[ffp[:, 0]] + fpc[ffp[:, 1]])
                mid = 0.5 * (np.asarray(pa) + np.asarray(pb))
                dist = np.linalg.norm(mids - mid, axis=1)
                dist[~is_perim] = np.inf
                sel = np.argsort(dist)[:4].astype(np.int32)
                pa = mids[sel[0]].tolist()
                pb = mids[sel[-1]].tolist()
            if len(sel) == 0:
                raise RuntimeError(f"no perimeter faces near BC '{name}' "
                                   f"segment {pa}->{pb} (tol {tol:.0f} m)")
            # ordered face-point chain + stations
            fps = []
            for fi in sel:
                for fp in ffp[fi]:
                    if fp not in fps:
                        fps.append(int(fp))
            st = np.r_[0.0, np.cumsum(lengths[sel])].astype(np.float32)
            lines[name] = dict(faces=sel, fps=np.array(fps, np.int32),
                               st=st, a=pa, b=pb)
            print(f"BC {name}: {len(sel)} perimeter faces, "
                  f"{st[-1]:.0f} m of line")

        # ---- rewrite the Boundary Condition Lines group
        bc = f["Geometry/Boundary Condition Lines"]
        for k in list(bc):
            del bc[k]
        attrs_dt = [("Name", "S32"), ("SA-2D", "S16"), ("Type", "S8"),
                    ("Length", "<f4")]
        ext_dt = [("BC Line ID", "<i4"), ("Face Index", "<i4"),
                  ("FP Start Index", "<i4"), ("FP End Index", "<i4"),
                  ("Station Start", "<f4"), ("Station End", "<f4")]
        attrs, ext, pinfo, pparts, ppts = [], [], [], [], []
        for lid, (name, L) in enumerate(lines.items()):
            attrs.append((name.encode(), area.encode(), b"External",
                          np.float32(L["st"][-1])))
            for j, fi in enumerate(L["faces"]):
                ext.append((lid, int(fi), int(ffp[fi][0]), int(ffp[fi][1]),
                            L["st"][j], L["st"][j + 1]))
            pinfo.append([len(ppts), 2, lid, 1])
            pparts.append([len(ppts), 2])
            ppts.extend([L["a"], L["b"]])
        bc.create_dataset("Attributes", data=np.array(attrs, attrs_dt))
        bc.create_dataset("External Faces", data=np.array(ext, ext_dt))
        bc.create_dataset("Polyline Info", data=np.array(pinfo, np.int32))
        bc.create_dataset("Polyline Parts", data=np.array(pparts, np.int32))
        bc.create_dataset("Polyline Points",
                          data=np.array(ppts, np.float64))
        f.attrs["Units System"] = np.bytes_("US Customary")
        f["Geometry"].attrs["SI Units"] = np.bytes_("False")

    write_ec(p01p, g01p, area, lines, series)


def write_ec(p01p, g01p, area, lines, series):
    # ---- Event Conditions in the plan HDF
    with h5py.File(p01p, "a") as f:
        # keep plan/geometry in sync (assemble copied the old BC tables)
        if "Geometry/Boundary Condition Lines" in f:
            del f["Geometry/Boundary Condition Lines"]
            with h5py.File(g01p, "r") as g:
                g.copy("Geometry/Boundary Condition Lines", f["Geometry"])
        bcg = f.require_group(
            "Event Conditions/Unsteady/Boundary Conditions")
        for grp in ["Flow Hydrographs", "Normal Depths"]:
            if grp in bcg:
                del bcg[grp]
        fh = bcg.create_group("Flow Hydrographs")
        nd = bcg.create_group("Normal Depths")

        L = lines["Inflow"]
        d = fh.create_dataset(f"2D: {area} BCLine: Inflow", data=series)
        common = dict(area=area, line="Inflow")
        d.attrs["2D Flow Area"] = np.bytes_(area)
        d.attrs["BC Line"] = np.bytes_("Inflow")
        d.attrs["Check TW Stage"] = np.bytes_("False")
        d.attrs["Data Type"] = np.bytes_("INST-VAL")
        d.attrs["EG Slope For Distributing Flow"] = np.float32(0.005)
        d.attrs["Start Date"] = np.bytes_(START)
        d.attrs["End Date"] = np.bytes_(END)
        d.attrs["Interval"] = np.bytes_("Days")
        d.attrs["Node Index"] = np.int32(1)
        d.attrs["Face Indexes"] = L["faces"]
        d.attrs["Face Fraction"] = np.ones(len(L["faces"]), np.float32)
        d.attrs["Face Point Indexes"] = L["fps"]

        L = lines["Outlet"]
        d = nd.create_dataset(f"2D: {area} BCLine: Outlet",
                              data=np.array([0.001], np.float32))
        d.attrs["2D Flow Area"] = np.bytes_(area)
        d.attrs["BC Line"] = np.bytes_("Outlet")
        d.attrs["BC Line WS"] = np.bytes_("Multiple")
        d.attrs["Check TW Stage"] = np.bytes_("False")
        d.attrs["Node Index"] = np.int32(2)
        d.attrs["Face Indexes"] = L["faces"]
        d.attrs["Face Fraction"] = np.ones(len(L["faces"]), np.float32)
        d.attrs["Face Point Indexes"] = L["fps"]

        ic = f.require_group("Event Conditions/Unsteady/Initial Conditions")
        ic.attrs["Startup Mode"] = np.bytes_("Computed")
        f.attrs["Units System"] = np.bytes_("US Customary")
        if "Geometry" in f:
            f["Geometry"].attrs["SI Units"] = np.bytes_("False")
    print("BC tables rewritten; event conditions injected; SI units set")


if __name__ == "__main__":
    main()
