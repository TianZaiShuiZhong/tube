from __future__ import annotations

import re
from pathlib import Path

from ..model import (
    Constraint,
    Element,
    ElementPressure,
    Material,
    Model,
    Node,
    NodalLoad,
    SectionPipe,
)


_NUM_RE = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")


def _parse_float(token: str) -> float:
    return float(token.strip())


def _split_csv_args(line: str) -> list[str]:
    # ANSYS command lines are CSV-like; keep empty fields.
    return [part.strip() for part in line.split(",")]


def parse_cdb(path: Path) -> Model:
    text = path.read_text(errors="ignore").splitlines()
    model = Model()

    # offsets (NUMOFF) exist in these files, but for MVP we keep IDs as-is.
    etype_map: dict[int, int] = {}  # type_num -> element_code (288/290)

    i = 0
    next_elem_id = 1

    while i < len(text):
        line = text[i].strip()
        if not line or line.startswith("/COM"):
            i += 1
            continue

        if line.startswith("ET,"):
            # ET, type_num, elem_code
            args = _split_csv_args(line)
            try:
                type_num = int(args[1])
                code = int(args[2])
                etype_map[type_num] = code
            except Exception:
                pass
            i += 1
            continue

        if line.startswith("NBLOCK"):
            # Skip format line then read node records until 'N,R5.3,LOC,-1,'
            i += 1
            if i < len(text) and text[i].strip().startswith("("):
                i += 1

            while i < len(text):
                l = text[i].rstrip("\n")
                if l.strip().startswith("N,R5.3"):
                    i += 1
                    break

                # Node record: id, ?, ?, x, y, z (y/z may be omitted)
                nums = _NUM_RE.findall(l)
                if len(nums) >= 4:
                    nid = int(nums[0])
                    x = float(nums[3])
                    y = float(nums[4]) if len(nums) >= 5 else 0.0
                    z = float(nums[5]) if len(nums) >= 6 else 0.0
                    model.nodes[nid] = Node(nid, x, y, z)

                i += 1
            continue

        if line.startswith("EBLOCK"):
            # After format line, each subsequent line describes one element in these benchmark CDBs.
            i += 1
            if i < len(text) and text[i].strip().startswith("("):
                i += 1

            while i < len(text):
                l = text[i].strip()
                if l == "-1":
                    i += 1
                    break

                ints = [int(x) for x in l.split() if x.strip()]
                if not ints:
                    i += 1
                    continue

                # Heuristic reverse-engineered from the provided benchmark CDBs:
                # - Common leading fields are attributes (mat/type/real/sec...).
                # - PIPE288 trailing fields are: EID, NI, NJ
                # - ELBOW290 trailing fields are: EID, NORIENT, NI, NJ
                type_num = ints[1] if len(ints) > 1 else 1
                elem_code = etype_map.get(type_num)

                if elem_code == 288 and len(ints) >= 3:
                    eid, n1, n2 = ints[-3], ints[-2], ints[-1]
                    model.elements[eid] = Element(
                        id=eid,
                        etype=288,
                        n1=n1,
                        n2=n2,
                        n3=None,
                        mat=ints[0] if len(ints) > 0 else None,
                        sec=ints[3] if len(ints) > 3 else None,
                    )
                elif elem_code == 290 and len(ints) >= 4:
                    # Tail fields: EID, NI, NJ, NORIENT (derived from benchmark files)
                    eid, n1, n2, norient = ints[-4], ints[-3], ints[-2], ints[-1]
                    model.elements[eid] = Element(
                        id=eid,
                        etype=290,
                        n1=n1,
                        n2=n2,
                        n3=norient,
                        mat=ints[0] if len(ints) > 0 else None,
                        sec=ints[3] if len(ints) > 3 else None,
                    )
                # else: unsupported element record

                i += 1
            continue

        if line.startswith("MPDATA"):
            # MPDATA,R5.0, 1,EX  , mat, 1, value
            args = _split_csv_args(line)
            if len(args) >= 7:
                prop = args[3].strip().upper()
                mat_id = int(args[4])
                value = _parse_float(args[6])
                old = model.materials.get(mat_id)
                if old is None:
                    old = Material(id=mat_id, ex=0.0, nuxy=0.0)

                if prop in {"EX", "EY", "EZ"}:
                    model.materials[mat_id] = Material(
                        id=mat_id,
                        ex=value,
                        nuxy=old.nuxy,
                        dens=old.dens,
                        alp=old.alp,
                        reft=old.reft,
                        creep=old.creep,
                    )
                elif prop in {"NUXY", "PRXY"}:
                    model.materials[mat_id] = Material(
                        id=mat_id,
                        ex=old.ex,
                        nuxy=value,
                        dens=old.dens,
                        alp=old.alp,
                        reft=old.reft,
                        creep=old.creep,
                    )
                elif prop in {"DENS"}:
                    model.materials[mat_id] = Material(
                        id=mat_id,
                        ex=old.ex,
                        nuxy=old.nuxy,
                        dens=value,
                        alp=old.alp,
                        reft=old.reft,
                        creep=old.creep,
                    )
                elif prop in {"ALPX"}:
                    model.materials[mat_id] = Material(
                        id=mat_id,
                        ex=old.ex,
                        nuxy=old.nuxy,
                        dens=old.dens,
                        alp=value,
                        reft=old.reft,
                        creep=old.creep,
                    )
                elif prop in {"REFT"}:
                    model.materials[mat_id] = Material(
                        id=mat_id,
                        ex=old.ex,
                        nuxy=old.nuxy,
                        dens=old.dens,
                        alp=old.alp,
                        reft=value,
                        creep=old.creep,
                    )
            i += 1
            continue

        if line.startswith("TB,BISO"):
            # TB,BISO, mat_id, ...
            args = _split_csv_args(line)
            mat_id = int(args[2]) if len(args) > 2 and args[2] else 1
            sy = etan = None

            j = i + 1
            while j < len(text):
                l2 = text[j].strip()
                if not l2:
                    j += 1
                    continue
                if l2.startswith("TB,") and j != i:
                    break
                if l2.startswith("TBDATA"):
                    nums = _NUM_RE.findall(l2)
                    vals = [float(x) for x in nums[1:]]
                    if vals:
                        if sy is None and len(vals) >= 1:
                            sy = vals[0]
                        if etan is None and len(vals) >= 2:
                            etan = vals[1]
                if l2.startswith("!END"):
                    break
                if l2.startswith("EXTOPT") or l2.startswith("TREF"):
                    break
                j += 1

            if sy is not None:
                if etan is None:
                    etan = 0.0
                old = model.materials.get(mat_id, Material(id=mat_id, ex=0.0, nuxy=0.0))
                model.materials[mat_id] = Material(
                    id=mat_id,
                    ex=old.ex,
                    nuxy=old.nuxy,
                    dens=old.dens,
                    alp=old.alp,
                    reft=old.reft,
                    creep=old.creep,
                    plasticity=(float(sy), float(etan)),
                )
            i = j
            continue

        if line.startswith("TB,CREE"):
            # Example:
            # TB,CREE,       1,   1,   4,1
            args = _split_csv_args(line)
            mat_id = int(args[2]) if len(args) > 2 and args[2] else 1
            c1 = c2 = c3 = c4 = None

            j = i + 1
            while j < len(text):
                l2 = text[j].strip()
                if not l2:
                    j += 1
                    continue
                if l2.startswith("TB,") and j != i:
                    break
                if l2.startswith("TBDATA"):
                    nums = _NUM_RE.findall(l2)
                    # nums[0] is the starting index (1 or 4)
                    vals = [float(x) for x in nums[1:]]
                    if vals:
                        if c1 is None and len(vals) >= 1:
                            c1 = vals[0]
                        if c2 is None and len(vals) >= 2:
                            c2 = vals[1]
                        if c3 is None and len(vals) >= 3:
                            c3 = vals[2]
                        if c4 is None and len(vals) >= 1 and nums and nums[0] == "4":
                            c4 = vals[0]
                if l2.startswith("!END"):
                    break
                # Stop when we reach something clearly unrelated
                if l2.startswith("EXTOPT") or l2.startswith("TREF"):
                    break
                j += 1

            if c1 is not None and c2 is not None and c3 is not None:
                if c4 is None:
                    c4 = 0.0
                old = model.materials.get(mat_id, Material(id=mat_id, ex=0.0, nuxy=0.0))
                model.materials[mat_id] = Material(
                    id=mat_id,
                    ex=old.ex,
                    nuxy=old.nuxy,
                    dens=old.dens,
                    alp=old.alp,
                    reft=old.reft,
                    creep=(float(c1), float(c2), float(c3), float(c4)),
                )

            i = j
            continue

        if line.startswith("NSUBST"):
            args = _split_csv_args(line)
            if len(args) >= 2 and args[1]:
                try:
                    model.nsubsteps = int(float(args[1]))
                except Exception:
                    pass
            i += 1
            continue

        if line.startswith("TIME"):
            args = _split_csv_args(line)
            if len(args) >= 2 and args[1]:
                try:
                    model.time_end = float(args[1])
                except Exception:
                    pass
            i += 1
            continue

        if line.startswith("SECTYPE"):
            # SECTYPE, sec_id, PIPE, , STRAI/BEND
            args = _split_csv_args(line)
            sec_id = int(args[1]) if len(args) > 1 and args[1] else 1
            # Next line is SECDATA
            i += 1
            if i < len(text) and text[i].strip().startswith("SECDATA"):
                sec_args = _split_csv_args(text[i].strip())
                # SECDATA, OD, t, ... (we only use first two)
                od = float(sec_args[1]) if len(sec_args) > 1 and sec_args[1] else 0.0
                t = float(sec_args[2]) if len(sec_args) > 2 and sec_args[2] else 0.0
                model.sections[sec_id] = SectionPipe(id=sec_id, outer_diameter=od, thickness=t)
            i += 1
            continue

        if line.startswith("BFUNIF,TEMP"):
            args = _split_csv_args(line)
            if len(args) >= 3:
                model.uniform_temperature = float(args[2])
            i += 1
            continue

        if line.startswith("ACEL"):
            args = _split_csv_args(line)
            if len(args) >= 4:
                ax = float(args[1]) if args[1] else 0.0
                ay = float(args[2]) if args[2] else 0.0
                az = float(args[3]) if args[3] else 0.0
                model.accel = (ax, ay, az)
            i += 1
            continue

        if line.startswith("D,"):
            # D, node, DOF, value
            args = _split_csv_args(line)
            if len(args) >= 4:
                node = int(args[1])
                dof = args[2].strip().upper()
                value = float(args[3]) if args[3] else 0.0
                model.constraints.append(Constraint(node=node, dof=dof, value=value))
            i += 1
            continue

        if line.startswith("F,"):
            # F, node, DOF, value
            args = _split_csv_args(line)
            if len(args) >= 4:
                node = int(args[1])
                dof = args[2].strip().upper()
                value = float(args[3]) if args[3] else 0.0
                model.nodal_loads.append(NodalLoad(node=node, dof=dof, value=value))
            i += 1
            continue

        if line.startswith("SFE,"):
            # SFE, elem, face, PRES, ... then next line contains magnitude
            args = _split_csv_args(line)
            if len(args) >= 2:
                try:
                    elem = int(args[1])
                    # pressure magnitude on next line
                    if i + 1 < len(text):
                        nums = _NUM_RE.findall(text[i + 1])
                        if nums:
                            p = float(nums[0])
                            model.elem_pressures.append(ElementPressure(elem=elem, value=p))
                            i += 2
                            continue
                except Exception:
                    pass
            i += 1
            continue

        i += 1

    return model
