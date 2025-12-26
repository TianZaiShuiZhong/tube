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
    i += 1  # LOOKUP_TABLE line
    i += 1
    # read node_ids
    for k in range(len(disp)):
        node_ids.append(int(lines[i + k].strip()))
    return node_ids, disp


def read_ansys_csv(path: Path):
    # Expect CSV: node_id, ux, uy, uz
    data = {}
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                nid = int(row[0])
                ux = float(row[1])
                uy = float(row[2])
                uz = float(row[3])
                data[nid] = (ux, uy, uz)
            except Exception:
                # skip header or malformed
                continue
    return data


def main():
    ap = argparse.ArgumentParser(description="Compare VTK displacement vs ANSYS CSV by node_id")
    ap.add_argument("--vtk", type=Path, required=True, help="VTK file exported by tubo (contains node_id scalar)")
    ap.add_argument("--ansys", type=Path, required=True, help="ANSYS CSV with columns: node_id,ux,uy,uz")
    args = ap.parse_args()

    nids, disp = read_vtk_points(args.vtk)
    ansys = read_ansys_csv(args.ansys)

    matched = 0
    err_sum = [0.0, 0.0, 0.0]
    err_max = [0.0, 0.0, 0.0]

    for nid, (ux, uy, uz) in zip(nids, disp):
        if nid in ansys:
            ax, ay, az = ansys[nid]
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
        print("No matching node_ids found.")
        return 1

    mae = [e / matched for e in err_sum]
    print(f"Compared {matched} nodes.")
    print(f"MAE ux,uy,uz: {mae[0]:.6g}, {mae[1]:.6g}, {mae[2]:.6g}")
    print(f"Max abs err ux,uy,uz: {err_max[0]:.6g}, {err_max[1]:.6g}, {err_max[2]:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
