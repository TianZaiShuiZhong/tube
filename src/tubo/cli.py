from __future__ import annotations

import argparse
from pathlib import Path

from .cdb.parser import parse_cdb
from .solver.creep import solve_creep_time_series
from .solver.linear_static import solve_linear_static
from .vtk.writer import write_vtk_polydata
from .vtk.pvd import write_pvd
from .vtk.surface import write_pipe_surface_quads


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tubo", description="Pipe-element solver (MVP)")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Parse CDB, run a linear static solve, export VTK")
    run.add_argument("--cdb", type=Path, required=True, help="ANSYS CDB input file")
    run.add_argument("--out", type=Path, required=True, help="VTK output path (.vtk)")
    run.add_argument("--surface", type=Path, default=None, help="Optional quad surface VTK output (.vtk)")
    run.add_argument("--circ", type=int, default=24, help="Circumferential divisions for surface quads")
    run.add_argument("--pressure-equiv", action="store_true", help="Enable pressure equivalence in assembly")

    creep = sub.add_parser("creep", help="Run a minimal creep time series and export a .pvd")
    creep.add_argument("--cdb", type=Path, required=True, help="ANSYS CDB input file")
    creep.add_argument("--outdir", type=Path, required=True, help="Output directory")
    creep.add_argument("--basename", type=str, default="result", help="Base name for outputs")
    creep.add_argument("--time", type=float, default=None, help="End time (override TIME in CDB)")
    creep.add_argument("--nsubsteps", type=int, default=None, help="Number of substeps (override NSUBST)")
    creep.add_argument("--pressure-equiv", action="store_true", help="Enable pressure equivalence in assembly")

    args = parser.parse_args(argv)

    if args.cmd == "run":
        model = parse_cdb(args.cdb)
        result = solve_linear_static(model, enable_pressure_equiv=bool(args.pressure_equiv))
        write_vtk_polydata(args.out, model, result)
        # optional surface expansion: use model section radius and centerline order from result node_ids
        if args.surface is not None:
            # try to fetch outer radius from model section
            # simplistic: SectionPipe has od and t (outer diameter and thickness)
            try:
                sec = next(iter(model.sections.values()))
                # SectionPipe stores `outer_diameter`; use its `ro` property when available
                radius = float(getattr(sec, "ro", 0.5 * getattr(sec, "outer_diameter", 1.0)))
            except Exception:
                radius = 1.0
            centerline = result["node_ids"].tolist()
            # Pass modal_state if available (currently None for linear_static)
            modal_state = model.modal_state
            write_pipe_surface_quads(args.surface, model, centerline, radius=radius, n_circ=int(args.circ), modal_state=modal_state)
        return 0

    if args.cmd == "creep":
        model = parse_cdb(args.cdb)
        outdir: Path = args.outdir
        outdir.mkdir(parents=True, exist_ok=True)

        time_end = args.time if args.time is not None else model.time_end
        nsubsteps = args.nsubsteps if args.nsubsteps is not None else model.nsubsteps
        if time_end is None or nsubsteps is None:
            raise SystemExit("Creep run requires TIME and NSUBST (either in CDB or via flags)")

        series = solve_creep_time_series(model, time_end=time_end, nsubsteps=nsubsteps, enable_pressure_equiv=bool(args.pressure_equiv))
        vtk_files: list[str] = []
        timesteps: list[float] = []
        for k, (t, res) in enumerate(series):
            vtk_name = f"{args.basename}_{k:04d}.vtk"
            vtk_path = outdir / vtk_name
            write_vtk_polydata(vtk_path, model, res)
            vtk_files.append(vtk_name)
            timesteps.append(float(t))

        write_pvd(outdir / f"{args.basename}.pvd", vtk_files=vtk_files, timesteps=timesteps)
        return 0

    raise RuntimeError("unreachable")
