from __future__ import annotations

from pathlib import Path


def write_pvd(path: Path, *, vtk_files: list[str], timesteps: list[float]) -> None:
    if len(vtk_files) != len(timesteps):
        raise ValueError("vtk_files and timesteps must have same length")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <Collection>\n")
        for t, name in zip(timesteps, vtk_files, strict=True):
            f.write(f'    <DataSet timestep="{t:.9g}" group="" part="0" file="{name}"/>\n')
        f.write("  </Collection>\n")
        f.write("</VTKFile>\n")
