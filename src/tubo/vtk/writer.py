from __future__ import annotations

from pathlib import Path

import numpy as np

from ..model import Model


def write_vtk_polydata(path: Path, model: Model, result: dict[str, np.ndarray]) -> None:
    node_ids = result["node_ids"].tolist()
    U = result["U"]

    # Build point list (in node_id order)
    pts = []
    disp = []
    for i, nid in enumerate(node_ids):
        n = model.nodes[nid]
        pts.append((n.x, n.y, n.z))
        disp.append((U[6 * i + 0], U[6 * i + 1], U[6 * i + 2]))

    # Build line cells from elements (n1-n2)
    lines = []
    id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}
    for e in model.elements.values():
        n1, n2 = e.end_nodes()
        if n1 in id_to_idx and n2 in id_to_idx:
            lines.append((id_to_idx[n1], id_to_idx[n2]))

    # Legacy VTK PolyData (ASCII)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("tubo result\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        f.write(f"POINTS {len(pts)} float\n")
        for x, y, z in pts:
            f.write(f"{x:.9g} {y:.9g} {z:.9g}\n")

        # LINES: each line is '2 i j'
        nlines = len(lines)
        f.write(f"LINES {nlines} {nlines * 3}\n")
        for i, j in lines:
            f.write(f"2 {i} {j}\n")

        f.write(f"POINT_DATA {len(pts)}\n")
        f.write("VECTORS displacement float\n")
        for ux, uy, uz in disp:
            f.write(f"{ux:.9g} {uy:.9g} {uz:.9g}\n")
        # also write node_id scalar for comparison scripts
        f.write("SCALARS node_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for nid in node_ids:
            f.write(f"{nid}\n")
