"""
Modal field expansion and visualization support for circumferential Fourier modes.
Scaffolding for m=0,1,2 (breathing, bending, oval) modes on pipe elements.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..model import Model, FourierModalState, Node


def compute_modal_displacement_at_angle(
    node_id: int,
    theta: float,
    modal_amplitudes: dict[int, float],
    radial: float = 1.0,
) -> tuple[float, float, float]:
    """
    Compute radial and circumferential displacement at a given angle theta
    on the pipe surface for a node with modal amplitude state.
    
    Args:
        node_id: Node identifier
        theta: Circumferential angle (0 to 2pi)
        modal_amplitudes: {mode: amplitude} dict for this node
        radial: Radial distance from pipe centerline
    
    Returns:
        (u_r, u_theta, u_z) displacement components in local cylindrical coords
    """
    u_r = 0.0
    u_t = 0.0
    u_z = 0.0
    
    for mode, amp in modal_amplitudes.items():
        if mode == 0:
            # m=0: breathing (radial breathing, axisymmetric)
            u_r += amp * radial  # simple radial scaling
        elif mode == 1:
            # m=1: bending-like (one lobe)
            u_r += amp * radial * math.cos(theta)
            u_t += amp * radial * math.sin(theta) * 0.5  # small tangential component
        elif mode == 2:
            # m=2: oval (two lobes)
            u_r += amp * radial * math.cos(2 * theta)
            u_t += amp * radial * math.sin(2 * theta) * 0.5
        # Additional modes can follow the pattern
    
    return u_r, u_t, u_z


def expand_surface_with_modal(
    model: Model,
    centerline_nodes: Sequence[int],
    radius: float,
    n_circ: int = 24,
    modal_state: FourierModalState | None = None,
) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int, int]], list[float]]:
    """
    Expand pipe centerline into 3D surface quads, optionally incorporating modal displacements.
    
    Args:
        model: Pipe model with nodes and modal_config
        centerline_nodes: Ordered list of centerline node IDs
        radius: Pipe outer radius
        n_circ: Circumferential discretization
        modal_state: Optional FourierModalState with modal amplitudes
    
    Returns:
        (points, quads, modal_radii)
        - points: 3D surface node coordinates (x, y, z)
        - quads: 4-node quad connectivity indices
        - modal_radii: Per-point radial displacement (for visualization coloring)
    """
    import numpy as np

    def _safe_dir(p0, p1):
        dx = np.array(p1) - np.array(p0)
        n = np.linalg.norm(dx)
        return (dx / n) if n > 0 else np.array([1.0, 0.0, 0.0])

    def _orthonormal(d):
        ax = np.array([0.0, 0.0, 1.0]) if abs(d[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(d, ax)
        nu = np.linalg.norm(u)
        u = (u / nu) if nu > 0 else np.array([1.0, 0.0, 0.0])
        v = np.cross(d, u)
        return u, v

    # Remove consecutive duplicate points to avoid zero-length tangents
    centers_raw = [model.nodes[nid] for nid in centerline_nodes]
    centers: list[Node] = []
    node_ids: list[int] = []
    tol = 1e-9
    for nid, nd in zip(centerline_nodes, centers_raw):
        if not centers:
            centers.append(nd)
            node_ids.append(nid)
            continue
        prev = centers[-1]
        if abs(prev.x - nd.x) < tol and abs(prev.y - nd.y) < tol and abs(prev.z - nd.z) < tol:
            continue
        centers.append(nd)
        node_ids.append(nid)

    # Corner smoothing: insert offset points to avoid pinching at sharp angles
    smoothed: list[Node] = []
    smoothed_ids: list[int] = []
    if centers:
        smoothed.append(centers[0])
        smoothed_ids.append(node_ids[0])
    for i in range(1, len(centers) - 1):
        c_prev, c_cur, c_next = centers[i - 1], centers[i], centers[i + 1]
        id_cur = node_ids[i]
        p_prev = np.array([c_prev.x, c_prev.y, c_prev.z])
        p_cur = np.array([c_cur.x, c_cur.y, c_cur.z])
        p_next = np.array([c_next.x, c_next.y, c_next.z])
        d_prev = _safe_dir(p_prev, p_cur)
        d_next = _safe_dir(p_cur, p_next)
        dot = float(np.clip(np.dot(d_prev, d_next), -1.0, 1.0))
        angle = math.acos(dot)
        len_prev = float(np.linalg.norm(p_cur - p_prev))
        len_next = float(np.linalg.norm(p_next - p_cur))
        smooth_len = min(0.5 * min(len_prev, len_next), radius)
        if angle > 1e-3 and smooth_len > 0:
            smoothed.append(Node(id=id_cur, x=c_cur.x - d_prev[0] * smooth_len, y=c_cur.y - d_prev[1] * smooth_len, z=c_cur.z - d_prev[2] * smooth_len))
            smoothed_ids.append(id_cur)
            smoothed.append(c_cur)
            smoothed_ids.append(id_cur)
            smoothed.append(Node(id=id_cur, x=c_cur.x + d_next[0] * smooth_len, y=c_cur.y + d_next[1] * smooth_len, z=c_cur.z + d_next[2] * smooth_len))
            smoothed_ids.append(id_cur)
        else:
            smoothed.append(c_cur)
            smoothed_ids.append(id_cur)
    if len(centers) > 1:
        smoothed.append(centers[-1])
        smoothed_ids.append(node_ids[-1])
    centers = smoothed
    node_ids = smoothed_ids

    if len(centers) < 2:
        return [], [], []

    # Build smooth, continuous frames along the path using parallel transport
    frames = []
    tangents: list[np.ndarray] = []
    for i in range(len(centers)):
        if i == 0:
            p0 = np.array([centers[i].x, centers[i].y, centers[i].z])
            p1 = np.array([centers[i + 1].x, centers[i + 1].y, centers[i + 1].z])
            d = _safe_dir(p0, p1)
        elif i == len(centers) - 1:
            p0 = np.array([centers[i - 1].x, centers[i - 1].y, centers[i - 1].z])
            p1 = np.array([centers[i].x, centers[i].y, centers[i].z])
            d = _safe_dir(p0, p1)
        else:
            p_prev = np.array([centers[i - 1].x, centers[i - 1].y, centers[i - 1].z])
            p_curr = np.array([centers[i].x, centers[i].y, centers[i].z])
            p_next = np.array([centers[i + 1].x, centers[i + 1].y, centers[i + 1].z])
            v_prev = p_curr - p_prev
            v_next = p_next - p_curr
            d_prev = _safe_dir(p_prev, p_curr)
            d_next = _safe_dir(p_curr, p_next)
            len_prev = float(np.linalg.norm(v_prev))
            len_next = float(np.linalg.norm(v_next))
            w_prev = len_prev if len_prev > 0 else 1.0
            w_next = len_next if len_next > 0 else 1.0
            d = w_prev * d_prev + w_next * d_next
            nd = np.linalg.norm(d)
            d = d_next if nd == 0 else (d / nd)
        tangents.append(d)

        if not frames:
            u, v = _orthonormal(d)
            frames.append((u, v))
        else:
            prev_u, prev_v = frames[-1]
            d_prev = tangents[-2]
            w = np.cross(d_prev, d)
            w_norm = float(np.linalg.norm(w))
            if w_norm < 1e-12:
                u, v = prev_u, prev_v
            else:
                w /= w_norm
                dot = float(np.clip(np.dot(d_prev, d), -1.0, 1.0))
                angle = math.acos(dot)

                def rotate(vec: np.ndarray) -> np.ndarray:
                    cos_a = math.cos(angle)
                    sin_a = math.sin(angle)
                    return (
                        vec * cos_a
                        + np.cross(w, vec) * sin_a
                        + w * (np.dot(w, vec)) * (1 - cos_a)
                    )

                u = rotate(prev_u)
                v = rotate(prev_v)
                u = u / np.linalg.norm(u)
                v = np.cross(d, u)
                v = v / np.linalg.norm(v)
            frames.append((u, v))

    pts = []
    modal_radii = []
    ring_offset = []
    ring_centers: list[tuple[float, float, float]] = []
    
    for i, c in enumerate(centers):
        u, v = frames[i]
        ring_offset.append(len(pts))
        ring_centers.append((c.x, c.y, c.z))
        
        nid = node_ids[i]
        modal_amps = {}
        if modal_state is not None:
            modal_amps = modal_state.modal_amplitudes.get(nid, {})
        
        for k in range(n_circ):
            theta = 2.0 * math.pi * (k / n_circ)
            
            # Compute modal displacement
            u_r, u_t, u_z = compute_modal_displacement_at_angle(
                nid, theta, modal_amps, radial=1.0
            )
            
            # Effective radius: base radius + modal radial displacement
            r_eff = radius * (1.0 + u_r / radius) if radius > 0 else radius
            
            # Surface point in local cylindrical coords
            px = c.x + r_eff * (math.cos(theta) * u[0] + math.sin(theta) * v[0])
            py = c.y + r_eff * (math.cos(theta) * u[1] + math.sin(theta) * v[1])
            pz = c.z + r_eff * (math.cos(theta) * u[2] + math.sin(theta) * v[2])
            
            pts.append((px, py, pz))
            modal_radii.append(float(u_r))  # Store radial modal component for coloring
    
    quads = []
    for i in range(len(centers) - 1):
        r0 = ring_offset[i]
        r1 = ring_offset[i + 1]
        c0 = np.array(ring_centers[i])
        c1 = np.array(ring_centers[i + 1])
        for k in range(n_circ):
            a = r0 + k
            b = r0 + ((k + 1) % n_circ)
            c_idx = r1 + ((k + 1) % n_circ)
            d_idx = r1 + k
            pa = np.array(pts[a])
            pb = np.array(pts[b])
            pc = np.array(pts[c_idx])
            pd = np.array(pts[d_idx])
            normal = np.cross(pb - pa, pc - pa)
            radial = (pa - c0) + (pb - c0) + (pc - c1) + (pd - c1)
            if float(np.dot(normal, radial)) < 0.0:
                quads.append((a, d_idx, c_idx, b))
            else:
                quads.append((a, b, c_idx, d_idx))

    return pts, quads, modal_radii
