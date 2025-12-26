from __future__ import annotations

import math

import numpy as np

from ..model import Model
from .assembly import assemble_linear_system


_DOF_ORDER = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")


def _section_props_pipe(od: float, t: float) -> tuple[float, float, float, float, float]:
    # Units follow CDB (mm, MPa, etc.).
    ro = 0.5 * od
    ri = ro - t
    if ri < 0:
        ri = 0.0
    area = math.pi * (ro * ro - ri * ri)
    iy = (math.pi / 4.0) * (ro**4 - ri**4)
    iz = iy
    j = (math.pi / 2.0) * (ro**4 - ri**4)  # polar (approx for circular tube)
    return ro, ri, area, iy, j


def _beam3d_local_k(E: float, G: float, A: float, Iy: float, Iz: float, J: float, L: float) -> np.ndarray:
    # 12x12 Euler-Bernoulli beam stiffness in local coordinates.
    k = np.zeros((12, 12), dtype=float)

    # axial
    k_ax = E * A / L
    k[0, 0] += k_ax
    k[0, 6] -= k_ax
    k[6, 0] -= k_ax
    k[6, 6] += k_ax

    # torsion (about local x)
    k_t = G * J / L
    k[3, 3] += k_t
    k[3, 9] -= k_t
    k[9, 3] -= k_t
    k[9, 9] += k_t

    # bending about local z (deflection in y, rotation about z)
    k1 = 12 * E * Iz / (L**3)
    k2 = 6 * E * Iz / (L**2)
    k3 = 4 * E * Iz / L
    k4 = 2 * E * Iz / L

    # dofs: y -> 1, rz -> 5
    k[1, 1] += k1
    k[1, 5] += k2
    k[1, 7] -= k1
    k[1, 11] += k2

    k[5, 1] += k2
    k[5, 5] += k3
    k[5, 7] -= k2
    k[5, 11] += k4

    k[7, 1] -= k1
    k[7, 5] -= k2
    k[7, 7] += k1
    k[7, 11] -= k2

    k[11, 1] += k2
    k[11, 5] += k4
    k[11, 7] -= k2
    k[11, 11] += k3

    # bending about local y (deflection in z, rotation about y)
    k1 = 12 * E * Iy / (L**3)
    k2 = 6 * E * Iy / (L**2)
    k3 = 4 * E * Iy / L
    k4 = 2 * E * Iy / L

    # dofs: z -> 2, ry -> 4
    k[2, 2] += k1
    k[2, 4] -= k2
    k[2, 8] -= k1
    k[2, 10] -= k2

    k[4, 2] -= k2
    k[4, 4] += k3
    k[4, 8] += k2
    k[4, 10] += k4

    k[8, 2] -= k1
    k[8, 4] += k2
    k[8, 8] += k1
    k[8, 10] += k2

    k[10, 2] -= k2
    k[10, 4] += k4
    k[10, 8] += k2
    k[10, 10] += k3

    return k


def _rotation_from_x_and_hint(x: np.ndarray, hint: np.ndarray) -> np.ndarray:
    # Build orthonormal basis: ex along x, ey from hint (projected), ez = ex x ey
    ex = x / np.linalg.norm(x)
    h = hint - np.dot(hint, ex) * ex
    if np.linalg.norm(h) < 1e-12:
        # pick a fallback
        tmp = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(tmp, ex)) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        h = tmp - np.dot(tmp, ex) * ex
    ey = h / np.linalg.norm(h)
    ez = np.cross(ex, ey)
    ez = ez / np.linalg.norm(ez)
    # local-to-global rotation matrix
    R = np.column_stack([ex, ey, ez])
    return R


def _beam_T(R: np.ndarray) -> np.ndarray:
    # 12x12 transform for a 3D beam with 6 dof per node
    T = np.zeros((12, 12), dtype=float)
    T[0:3, 0:3] = R
    T[3:6, 3:6] = R
    T[6:9, 6:9] = R
    T[9:12, 9:12] = R
    return T


def _udl_equiv_nodal_local(wy: float, wz: float, L: float) -> np.ndarray:
    # Uniform distributed load in local coordinates (forces per length) applied to y and z.
    # Returns 12 vector in local DOF order.
    f = np.zeros(12, dtype=float)

    # in local y
    f[1] += wy * L / 2
    f[7] += wy * L / 2
    f[5] += wy * (L**2) / 12
    f[11] -= wy * (L**2) / 12

    # in local z
    f[2] += wz * L / 2
    f[8] += wz * L / 2
    f[4] -= wz * (L**2) / 12
    f[10] += wz * (L**2) / 12

    return f


def solve_linear_static(model: Model, enable_pressure_equiv: bool = True) -> dict[str, np.ndarray]:
    # Assemble via shared assembly to keep behavior consistent across run/creep
    node_ids, node_index, K, F, _elem_infos, dof_map, _props = assemble_linear_system(
        model, extra_initial_strain_by_elem=None, enable_pressure_equiv=enable_pressure_equiv
    )

    ndof = K.shape[0]
    # Detect DOF count per node from K size
    dof_per_node = ndof // len(node_ids)
    
    fixed: dict[int, float] = {}
    for bc in model.constraints:
        if bc.node not in node_index:
            continue
        if bc.dof not in dof_map:
            continue
        fixed[dof_per_node * node_index[bc.node] + dof_map[bc.dof]] = bc.value

    free_dofs = np.array([i for i in range(ndof) if i not in fixed], dtype=int)
    fixed_dofs = np.array(sorted(fixed.keys()), dtype=int)

    U = np.zeros(ndof, dtype=float)
    if fixed_dofs.size:
        U[fixed_dofs] = np.array([fixed[i] for i in fixed_dofs], dtype=float)
    if free_dofs.size:
        Kff = K[np.ix_(free_dofs, free_dofs)]
        Ff = F[free_dofs].copy()
        if fixed_dofs.size:
            Kfc = K[np.ix_(free_dofs, fixed_dofs)]
            Ff -= Kfc @ U[fixed_dofs]
        Uf = np.linalg.solve(Kff, Ff)
        U[free_dofs] = Uf

    # Extract modal state if available (DOF 6 is OVAL)
    if dof_per_node >= 7 and "OVAL" in dof_map:
        oval_idx = dof_map["OVAL"]
        modal_amplitudes = {}
        for i, nid in enumerate(node_ids):
            amp = U[dof_per_node * i + oval_idx]
            # Store as mode 2 (ovalization)
            modal_amplitudes[nid] = {2: float(amp)}
        
        # Update model with new modal state (hacky side-effect, but effective for MVP)
        # Ideally we should return this in the result dict
        from ..model import FourierModalState
        model.modal_state = FourierModalState(modal_amplitudes=modal_amplitudes)

    return {"node_ids": np.array(node_ids, dtype=int), "U": U}
