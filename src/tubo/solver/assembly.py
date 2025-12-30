from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..model import Model, Material


_DOF_ORDER = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ", "OVAL")


@dataclass(frozen=True)
class ElementKinematics:
    eid: int
    n1: int
    n2: int
    L: float
    R: np.ndarray  # 3x3 local-to-global rotation
    T: np.ndarray  # 14x14
    dofs: np.ndarray  # 14 global dof indices


def _section_props_pipe(od: float, t: float) -> tuple[float, float, float, float, float]:
    ro = 0.5 * od
    ri = ro - t
    if ri < 0:
        ri = 0.0
    area = math.pi * (ro * ro - ri * ri)
    iy = (math.pi / 4.0) * (ro**4 - ri**4)
    iz = iy
    j = (math.pi / 2.0) * (ro**4 - ri**4)
    return ro, ri, area, iy, j


def _beam3d_local_k(E: float, G: float, A: float, Iy: float, Iz: float, J: float, L: float, 
                    ro: float, ri: float, nu: float, curvature: float = 0.0) -> np.ndarray:
    k = np.zeros((14, 14), dtype=float)

    # axial
    k_ax = E * A / L
    k[0, 0] += k_ax
    k[0, 7] -= k_ax
    k[7, 0] -= k_ax
    k[7, 7] += k_ax

    # torsion
    k_t = G * J / L
    k[3, 3] += k_t
    k[3, 10] -= k_t
    k[10, 3] -= k_t
    k[10, 10] += k_t

    # bending about local z (y deflection)
    k1 = 12 * E * Iz / (L**3)
    k2 = 6 * E * Iz / (L**2)
    k3 = 4 * E * Iz / L
    k4 = 2 * E * Iz / L
    k[1, 1] += k1
    k[1, 5] += k2
    k[1, 8] -= k1
    k[1, 12] += k2
    k[5, 1] += k2
    k[5, 5] += k3
    k[5, 8] -= k2
    k[5, 12] += k4
    k[8, 1] -= k1
    k[8, 5] -= k2
    k[8, 8] += k1
    k[8, 12] -= k2
    k[12, 1] += k2
    k[12, 5] += k4
    k[12, 8] -= k2
    k[12, 12] += k3

    # bending about local y (z deflection)
    k1 = 12 * E * Iy / (L**3)
    k2 = 6 * E * Iy / (L**2)
    k3 = 4 * E * Iy / L
    k4 = 2 * E * Iy / L
    k[2, 2] += k1
    k[2, 4] -= k2
    k[2, 9] -= k1
    k[2, 11] -= k2
    k[4, 2] -= k2
    k[4, 4] += k3
    k[4, 9] += k2
    k[4, 11] += k4
    k[9, 2] -= k1
    k[9, 4] += k2
    k[9, 9] += k1
    k[9, 11] += k2
    k[11, 2] -= k2
    k[11, 4] += k4
    k[11, 9] += k2
    k[11, 11] += k3

    # Ovalization (m=2) stiffness
    # Ring bending stiffness D = E * t^3 / (12 * (1 - nu^2))
    # K_oval = (9 * pi * D / R^3) * L
    t = ro - ri
    R_mean = (ro + ri) / 2.0
    D = E * (t**3) / (12.0 * (1.0 - nu**2))
    k_oval = (9.0 * math.pi * D / (R_mean**3)) * L
    
    if k_oval < 1e-6:
        k_oval = 1e-6

    k[6, 6] += k_oval / 2.0
    k[13, 13] += k_oval / 2.0
    
    # Karman effect coupling: Bending (about Z) <-> Ovalization
    # Only active if curvature > 0 (Elbow)
    # Coupling term C = 3 / (2 * R_curv)
    # Stiffness contribution: - E * I * C * (theta_z2 - theta_z1) * (oval1 + oval2)/2
    # This is a simplified interaction term.
    if abs(curvature) > 1e-9:
        C = 3.0 * curvature / 2.0
        # Coupling coefficient
        # K_c = E * Iz * C / 2 (factor 1/2 from averaging ovalization)
        # But wait, energy is - E I C kappa a.
        # kappa ~ (theta_z2 - theta_z1) / L
        # a ~ (a1 + a2) / 2
        # Energy U_c = - E * Iz * C * (th_z2 - th_z1)/L * (a1 + a2)/2 * L
        #            = - E * Iz * C * 0.5 * (th_z2 - th_z1) * (a1 + a2)
        # Force F_i = dU_c / d_i
        # F_th_z1 = + 0.5 * E * Iz * C * (a1 + a2)
        # F_th_z2 = - 0.5 * E * Iz * C * (a1 + a2)
        # F_a1    = - 0.5 * E * Iz * C * (th_z2 - th_z1)
        # F_a2    = - 0.5 * E * Iz * C * (th_z2 - th_z1)
        
        Kc = 0.5 * E * Iz * C
        
        # Rows: 5 (th_z1), 12 (th_z2), 6 (a1), 13 (a2)
        # F_th_z1 terms
        k[5, 6] += Kc
        k[5, 13] += Kc
        
        # F_th_z2 terms
        k[12, 6] -= Kc
        k[12, 13] -= Kc
        
        # F_a1 terms
        k[6, 5] += Kc
        k[6, 12] -= Kc
        
        # F_a2 terms
        k[13, 5] += Kc
        k[13, 12] -= Kc

    return k


def _interp_prop(mat: Material, key: str, temp: float | None, default: float | None):
    """Linear interpolation of temperature-dependent property."""
    if temp is None:
        return default
    table = mat.temp_tables.get(key)
    if not table:
        return default
    pts = sorted(table, key=lambda x: x[0])
    for t, v in pts:
        if abs(t - temp) < 1e-9:
            return v
    if temp <= pts[0][0]:
        return pts[0][1]
    if temp >= pts[-1][0]:
        return pts[-1][1]
    for (t0, v0), (t1, v1) in zip(pts[:-1], pts[1:]):
        if t0 <= temp <= t1:
            alpha = (temp - t0) / (t1 - t0)
            return v0 + alpha * (v1 - v0)
    return default


def _circle_from_three(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> tuple[float, np.ndarray, np.ndarray] | None:
    """Compute circle radius, center, and normal from three non-collinear points."""
    v1 = p2 - p1
    v2 = p3 - p1
    n = np.cross(v1, v2)
    n_norm = float(np.linalg.norm(n))
    if n_norm < 1e-9:
        return None
    n = n / n_norm
    v1_sq = np.dot(v1, v1)
    v2_sq = np.dot(v2, v2)
    # perpendicular bisector intersection in plane
    mat3 = np.column_stack([v1, v2, n])
    rhs = 0.5 * np.array([v1_sq, v2_sq, np.dot(n, p1)], dtype=float)
    try:
        sol = np.linalg.solve(mat3, rhs)
    except np.linalg.LinAlgError:
        return None
    center = p1 + sol[0] * v1 + sol[1] * v2
    radius = float(np.linalg.norm(p1 - center))
    if radius < 1e-9:
        return None
    return radius, center, n


def _rotation_from_x_and_hint(x: np.ndarray, hint: np.ndarray) -> np.ndarray:
    ex = x / np.linalg.norm(x)
    h = hint - np.dot(hint, ex) * ex
    if np.linalg.norm(h) < 1e-12:
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(tmp, ex)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        h = tmp - np.dot(tmp, ex) * ex
    ey = h / np.linalg.norm(h)
    ez = np.cross(ex, ey)
    ez = ez / np.linalg.norm(ez)
    return np.column_stack([ex, ey, ez])


def _beam_T(R: np.ndarray) -> np.ndarray:
    T = np.zeros((14, 14), dtype=float)
    T[0:3, 0:3] = R
    T[3:6, 3:6] = R
    # DOF 6 (OVAL) is scalar, identity
    T[6, 6] = 1.0
    
    T[7:10, 7:10] = R
    T[10:13, 10:13] = R
    # DOF 13 (OVAL) is scalar, identity
    T[13, 13] = 1.0
    return T


def _udl_equiv_nodal_local(wx: float, wy: float, wz: float, L: float) -> np.ndarray:
    f = np.zeros(14, dtype=float)

    f[0] += wx * L / 2
    f[7] += wx * L / 2

    f[1] += wy * L / 2
    f[8] += wy * L / 2
    f[5] += wy * (L**2) / 12
    f[12] -= wy * (L**2) / 12

    f[2] += wz * L / 2
    f[9] += wz * L / 2
    f[4] -= wz * (L**2) / 12
    f[11] += wz * (L**2) / 12

    return f


def collect_active_nodes(model: Model) -> list[int]:
    active_nodes: set[int] = set()
    for e in model.elements.values():
        n1, n2 = e.end_nodes()
        active_nodes.add(n1)
        active_nodes.add(n2)
    for bc in model.constraints:
        active_nodes.add(bc.node)
    for load in model.nodal_loads:
        active_nodes.add(load.node)
    return sorted(nid for nid in active_nodes if nid in model.nodes)


def assemble_linear_system(
    model: Model,
    *,
    extra_initial_strain_by_elem: dict[int, float] | None = None,
    extra_initial_curvature_by_elem: dict[int, tuple[float, float]] | None = None,
    enable_pressure_equiv: bool = True,
) -> tuple[list[int], dict[int, int], np.ndarray, np.ndarray, list[ElementKinematics], dict[str, int], tuple[float, float, float]]:
    node_ids = collect_active_nodes(model)
    node_index = {nid: idx for idx, nid in enumerate(node_ids)}

    ndof = 7 * len(node_ids)
    K = np.zeros((ndof, ndof), dtype=float)
    F = np.zeros(ndof, dtype=float)

    dof_map = {name: j for j, name in enumerate(_DOF_ORDER)}

    if not model.materials:
        raise ValueError("No materials parsed from CDB")
    mat = model.materials[min(model.materials.keys())]
    temp_use = model.uniform_temperature if model.uniform_temperature is not None else mat.reft
    ex_eff = _interp_prop(mat, "EX", temp_use, mat.ex)
    nu_eff = _interp_prop(mat, "NUXY", temp_use, mat.nuxy)
    alp_eff = _interp_prop(mat, "ALPX", temp_use, mat.alp)
    dens_eff = _interp_prop(mat, "DENS", temp_use, mat.dens)

    if not model.sections:
        raise ValueError("No sections parsed from CDB")
    sec = model.sections[min(model.sections.keys())]

    ro, ri, A, Iy, J = _section_props_pipe(sec.outer_diameter, sec.thickness)
    Iz = Iy

    E = ex_eff
    nu = nu_eff
    G = E / (2.0 * (1.0 + nu))
    ai = math.pi * (ri**2)

    p_by_elem: dict[int, float] = {}
    p_default: float | None = None
    if model.elem_pressures:
        p_by_elem = {ep.elem: ep.value for ep in model.elem_pressures}
        p_default = float(sum(ep.value for ep in model.elem_pressures) / len(model.elem_pressures))

    elem_infos: list[ElementKinematics] = []

    for e in model.elements.values():
        n1, n2 = e.end_nodes()
        if n1 not in node_index or n2 not in node_index:
            continue

        p1 = np.array([model.nodes[n1].x, model.nodes[n1].y, model.nodes[n1].z], dtype=float)
        p2 = np.array([model.nodes[n2].x, model.nodes[n2].y, model.nodes[n2].z], dtype=float)
        dx = p2 - p1
        L = float(np.linalg.norm(dx))
        if L <= 0:
            continue

        on = e.orientation_node()
        curvature = 0.0
        R = None

        if e.etype == 290 and on is not None and on in model.nodes:
            p3 = np.array([model.nodes[on].x, model.nodes[on].y, model.nodes[on].z], dtype=float)
            circle = _circle_from_three(p1, p2, p3)
            if circle is not None:
                r, center, normal = circle
                t1 = np.cross(normal, (p1 - center))
                t2 = np.cross(normal, (p2 - center))
                if np.linalg.norm(t1) > 1e-9 and np.linalg.norm(t2) > 1e-9:
                    t1 = t1 / np.linalg.norm(t1)
                    t2 = t2 / np.linalg.norm(t2)
                    dot = float(np.clip(np.dot(t1, t2), -1.0, 1.0))
                    ang = float(math.acos(dot))
                    if ang > 1e-6:
                        L = r * ang
                        curvature = 1.0 / r
                        t_mid = t1 + t2
                        if np.linalg.norm(t_mid) < 1e-9:
                            t_mid = t1
                        t_mid = t_mid / np.linalg.norm(t_mid)
                        b_vec = np.cross(normal, t_mid)
                        if np.linalg.norm(b_vec) < 1e-9:
                            b_vec = t_mid
                        b_vec = b_vec / np.linalg.norm(b_vec)
                        R = np.column_stack([t_mid, b_vec, normal])
        if R is None:
            hint = np.array([0.0, 0.0, 1.0])
            if on is not None and on in model.nodes:
                op = np.array([model.nodes[on].x, model.nodes[on].y, model.nodes[on].z], dtype=float)
                hint = op - p1
            R = _rotation_from_x_and_hint(dx, hint)

        T = _beam_T(R)

        k_local = _beam3d_local_k(E, G, A, Iy, Iz, J, L, ro, ri, nu, curvature=curvature)
        k_global = T.T @ k_local @ T

        dofs = []
        for nid in (n1, n2):
            base = 7 * node_index[nid]
            dofs.extend([base + j for j in range(7)])
        dofs_arr = np.array(dofs, dtype=int)

        for a in range(14):
            ia = dofs_arr[a]
            for b in range(14):
                ib = dofs_arr[b]
                K[ia, ib] += k_global[a, b]

        elem_infos.append(ElementKinematics(eid=e.id, n1=n1, n2=n2, L=L, R=R, T=T, dofs=dofs_arr))

        # self weight
        if model.accel is not None and dens_eff is not None:
            ax, ay, az = model.accel
            gvec = np.array([ax, ay, az], dtype=float)
            qg = dens_eff * A * gvec
            ql = R.T @ qg
            f_local = _udl_equiv_nodal_local(wx=ql[0], wy=ql[1], wz=ql[2], L=L)
            f_global = T.T @ f_local
            F[dofs_arr] += f_global

        # curved pressure load (Bourdon-like) per element when curvature is known
        if enable_pressure_equiv and curvature > 0.0 and model.elem_pressures:
            p_elem = None
            if e.id in p_by_elem:
                p_elem = p_by_elem[e.id]
            elif p_default is not None:
                p_elem = p_default
            if p_elem is not None:
                scale = mat.k_bourdon_scale if mat.k_bourdon_scale is not None else 1.0
                q_global = (-p_elem * ai * curvature * scale) * R[:, 2]  # normal points to center, force outward
                q_local = R.T @ q_global
                f_local = _udl_equiv_nodal_local(wx=q_local[0], wy=q_local[1], wz=q_local[2], L=L)
                f_global = T.T @ f_local
                F[dofs_arr] += f_global

        # initial strain: thermal + extra (creep)
        eps0 = 0.0
        if model.uniform_temperature is not None and alp_eff is not None and mat.reft is not None:
            dT = model.uniform_temperature - mat.reft
            eps0 += alp_eff * dT
        if extra_initial_strain_by_elem is not None:
            eps0 += float(extra_initial_strain_by_elem.get(e.id, 0.0))

        ky0, kz0 = 0.0, 0.0
        if extra_initial_curvature_by_elem is not None:
            ky0, kz0 = extra_initial_curvature_by_elem.get(e.id, (0.0, 0.0))

        if abs(eps0) > 0 or abs(ky0) > 0 or abs(kz0) > 0:
            f0_local = np.zeros(14, dtype=float)
            
            if abs(eps0) > 0:
                N0 = E * A * eps0
                f0_local[0] -= N0
                f0_local[7] += N0
            
            if abs(ky0) > 0:
                My0 = E * Iy * ky0
                # Bending about local y (affects rotation about y)
                # Equivalent nodal moments: -My0 at node 1, +My0 at node 2
                f0_local[4] -= My0
                f0_local[11] += My0

            if abs(kz0) > 0:
                Mz0 = E * Iz * kz0
                # Bending about local z (affects rotation about z)
                # Equivalent nodal moments: -Mz0 at node 1, +Mz0 at node 2
                f0_local[5] -= Mz0
                f0_local[12] += Mz0

            f0_global = T.T @ f0_local
            F[dofs_arr] += f0_global

    # nodal loads
    for load in model.nodal_loads:
        if load.node not in node_index:
            continue
        if load.dof not in dof_map:
            continue
        F[7 * node_index[load.node] + dof_map[load.dof]] += load.value

    # pressure handling: end-cap thrust (default) and future wall-pressure equivalence
    if enable_pressure_equiv and model.elem_pressures:
        adj: dict[int, list[int]] = {nid: [] for nid in node_ids}
        for e in model.elements.values():
            n1, n2 = e.end_nodes()
            if n1 in adj and n2 in adj:
                adj[n1].append(n2)
                adj[n2].append(n1)

        p_by_elem = {ep.elem: ep.value for ep in model.elem_pressures}
        ai = math.pi * (ri**2)

        for nid in node_ids:
            neigh = adj.get(nid, [])
            if len(neigh) != 1:
                continue
            other = neigh[0]

            pvals = []
            for e in model.elements.values():
                if {e.n1, e.n2} == {nid, other} and e.id in p_by_elem:
                    pvals.append(p_by_elem[e.id])
            if not pvals:
                pvals = [float(ep.value) for ep in model.elem_pressures]
            p = float(sum(pvals) / len(pvals))

            npos = np.array([model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z], dtype=float)
            opos = np.array([model.nodes[other].x, model.nodes[other].y, model.nodes[other].z], dtype=float)
            t = npos - opos
            tn = float(np.linalg.norm(t))
            if tn <= 0:
                continue
            t = t / tn

            fvec = p * ai * t
            base = 7 * node_index[nid]
            F[base + 0] += fvec[0]
            F[base + 1] += fvec[1]
            F[base + 2] += fvec[2]

        # Distributed load due to curvature (replacing corner resultant)
        # This converts internal pressure acting on curved pipe into equivalent
        # distributed loads on the beam elements.
        
        # Map edges to elements for pressure lookup
        edge_to_elem = {}
        for e in model.elements.values():
            edge = tuple(sorted((e.n1, e.n2)))
            edge_to_elem[edge] = e.id

        # Default pressure if specific element pressure is missing
        p_default = 0.0
        if model.elem_pressures:
            p_default = float(sum(ep.value for ep in model.elem_pressures) / len(model.elem_pressures))

        # Calculate distributed load vector q at each internal node
        node_q = {}  # nid -> np.array([qx, qy, qz])
        
        for nid in node_ids:
            neigh = adj.get(nid, [])
            if len(neigh) != 2:
                continue

            n1_id, n2_id = neigh
            p0 = np.array([model.nodes[nid].x, model.nodes[nid].y, model.nodes[nid].z], dtype=float)
            p1 = np.array([model.nodes[n1_id].x, model.nodes[n1_id].y, model.nodes[n1_id].z], dtype=float)
            p2 = np.array([model.nodes[n2_id].x, model.nodes[n2_id].y, model.nodes[n2_id].z], dtype=float)

            v1 = p1 - p0
            v2 = p2 - p0
            l1 = float(np.linalg.norm(v1))
            l2 = float(np.linalg.norm(v2))

            if l1 <= 1e-9 or l2 <= 1e-9:
                continue

            t1 = v1 / l1
            t2 = v2 / l2

            # If nearly colinear, curvature is zero
            if abs(np.dot(t1, t2)) > 0.9999:
                continue

            # Pressure magnitude
            e1_id = edge_to_elem.get(tuple(sorted((nid, n1_id))))
            e2_id = edge_to_elem.get(tuple(sorted((nid, n2_id))))
            
            val1 = p_by_elem.get(e1_id, p_default) if e1_id is not None else p_default
            val2 = p_by_elem.get(e2_id, p_default) if e2_id is not None else p_default
            p_local = (val1 + val2) / 2.0

            # Force direction and magnitude
            # Net force at node F = P * Ai * (-t1 - t2)
            # Distributed load q = F / L_tributary
            # L_tributary = (l1 + l2) / 2
            
            f_node = p_local * ai * (-t1 - t2)
            node_q[nid] = f_node / (0.5 * (l1 + l2))

        # Apply distributed loads to elements
        for info in elem_infos:
            q1 = node_curvature_force_val = node_q.get(info.n1, np.zeros(3))
            q2 = node_curvature_force_val = node_q.get(info.n2, np.zeros(3))

            if np.linalg.norm(q1) < 1e-12 and np.linalg.norm(q2) < 1e-12:
                continue

            # Average distributed load on element
            q_global = 0.5 * (q1 + q2)
            
            # Transform to local coordinates
            # info.R is 3x3 rotation matrix (columns are local axes in global coords)
            # q_local = R^T * q_global
            q_local = info.R.T @ q_global
            
            # q_local[0] is axial (x), q_local[1] is y, q_local[2] is z
            # We apply transverse loads and axial loads
            f_equiv_local = _udl_equiv_nodal_local(wx=q_local[0], wy=q_local[1], wz=q_local[2], L=info.L)
            
            # Transform forces back to global
            f_equiv_global = info.T.T @ f_equiv_local
            
            # Add to global force vector
            F[info.dofs] += f_equiv_global

    return node_ids, node_index, K, F, elem_infos, dof_map, (E, A, float(model.uniform_temperature or 0.0))
