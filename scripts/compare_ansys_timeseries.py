#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_vtk_points(path: Path):
    node_ids = []
    disp = []
    with path.open("r", encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f]
    # find POINT_DATA
    i = 0
    while i < len(lines) and not lines[i].startswith("POINT_DATA"):
        i += 1
    if i >= len(lines):
        raise ValueError("POINT_DATA not found in VTK")
    i += 1
    if i >= len(lines) or not lines[i].startswith("VECTORS displacement"):
        raise ValueError("displacement vectors not found in VTK")
    i += 1
    # read displacement vectors until SCALARS
    while i < len(lines) and not lines[i].startswith("SCALARS"):
        parts = lines[i].split()
        if len(parts) == 3:
            ux, uy, uz = map(float, parts)
            disp.append((ux, uy, uz))
        i += 1
    if i >= len(lines) or not lines[i].startswith("SCALARS node_id"):
        raise ValueError("node_id scalars not found in VTK (update writer)")
    i += 1  # LOOKUP_TABLE default
    i += 1
    # read node_ids
    for k in range(len(disp)):
        node_ids.append(int(lines[i + k].strip()))
    return node_ids, disp


def read_pvd(path: Path):
    # Very simple PVD reader: parse XML lines with file and timestep
    vtk_files = []
    timesteps = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith("<DataSet") and "file=" in line and "timestep=" in line:
            # extract attributes
            def get_attr(s: str, name: str) -> str:
                i = s.find(name + "=\"")
                if i == -1:
                    return ""
                j = s.find("\"", i + len(name) + 2)
                return s[i + len(name) + 2 : j]
            file = get_attr(line, "file")
            ts = get_attr(line, "timestep")
            vtk_files.append(file)
            timesteps.append(float(ts))
    return timesteps, vtk_files


def read_ansys_csv_series(dir_path: Path):
    # Expect files named like step_0000.csv each with columns: time,node_id,ux,uy,uz
    series = []
    for csv_path in sorted(dir_path.glob("*.csv")):
        time_val = None
        data = {}
        with csv_path.open("r", encoding="utf-8") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                try:
                    # time can be in first column or a header row; prefer a 'time' column
                    if row[0].lower() == "time":
                        time_val = float(row[1])
                        continue
                    nid = int(row[0])
                    ux = float(row[1])
                    uy = float(row[2])
                    uz = float(row[3])
                    data[nid] = (ux, uy, uz)
                except Exception:
                    continue
        # derive time from filename if not found
        if time_val is None:
            try:
                # look for pattern like _123.45 in name
                s = csv_path.stem
                for part in s.split("_"):
                    try:
                        time_val = float(part)
                        break
                    except Exception:
                        pass
            except Exception:
                time_val = 0.0
        series.append((float(time_val or 0.0), data, csv_path.name))
    return series


def main():
    ap = argparse.ArgumentParser(description="Compare VTK .pvd time series vs ANSYS CSV series")
    ap.add_argument("--pvd", type=Path, required=True, help="PVD file indexed to VTK timesteps")
    ap.add_argument("--vtkdir", type=Path, required=True, help="Directory containing VTK files referenced by PVD")
    ap.add_argument("--ansysdir", type=Path, required=True, help="Directory containing ANSYS CSV files (time,node_id,ux,uy,uz)")
    args = ap.parse_args()

    ts_pvd, files_pvd = read_pvd(args.pvd)
    ansys_series = read_ansys_csv_series(args.ansysdir)

    # build mapping from time to ANSYS data (nearest match)
    ansys_series.sort(key=lambda x: x[0])

    def nearest_ansys(time_val: float):
        if not ansys_series:
            return None
        best = None
        best_diff = 1e99
        for t, data, name in ansys_series:
            d = abs(t - time_val)
            if d < best_diff:
                best_diff = d
                best = (t, data, name)
        return best

    print(f"PVD steps: {len(ts_pvd)}; ANSYS steps: {len(ansys_series)}")
    for t, vtk_name in zip(ts_pvd, files_pvd):
        vtk_path = args.vtkdir / vtk_name
        nids, disp = read_vtk_points(vtk_path)
        ansys = nearest_ansys(t)
        if ansys is None:
            print(f"time {t:.6g}: no ANSYS data")
            continue
        t_ansys, data_ansys, name = ansys
        matched = 0
        err_sum = [0.0, 0.0, 0.0]
        err_max = [0.0, 0.0, 0.0]
        for nid, (ux, uy, uz) in zip(nids, disp):
            if nid in data_ansys:
                ax, ay, az = data_ansys[nid]
                ex = ux - ax
                ey = uy - ay
                ez = uz - az
                err_sum[0] += abs(ex)
                err_sum[1] += abs(ey)
                err_sum[2] += abs(ez)
                err_max[0] = max(err_max[0], abs(ex))
                err_max[1] = max(err_max[1], abs(ey))
                err_max[2] = max(err_max[2], abs(ez))
                matched += 1
        if matched == 0:
            print(f"time {t:.6g}: matched 0 nodes (ANSYS file {name})")
            continue
        mae = [e / matched for e in err_sum]
        print(f"time {t:.6g} (ANSYS {t_ansys:.6g} {name}): matched {matched}, MAE ux,uy,uz = {mae[0]:.6g}, {mae[1]:.6g}, {mae[2]:.6g}; max = {err_max[0]:.6g}, {err_max[1]:.6g}, {err_max[2]:.6g}")


if __name__ == "__main__":
    raise SystemExit(main())
