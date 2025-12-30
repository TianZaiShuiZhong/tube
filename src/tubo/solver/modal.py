"""
Modal field expansion and visualization support for circumferential Fourier modes.
Scaffolding for m=0,1,2 (breathing, bending, oval) modes on pipe elements.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..model import Model, FourierModalState


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
    centers: list = []
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

    if len(centers) < 2:
        return [], [], []

    # Build smooth, continuous frames along the path to avoid orientation flips
    frames = []
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

        u, v = _orthonormal(d)
        if frames:
            prev_u, _ = frames[-1]
            if float(np.dot(prev_u, u)) < 0.0:
                u = -u
                v = -v
        frames.append((u, v))

    pts = []
    modal_radii = []
    ring_offset = []
    
    for i, c in enumerate(centers):
        u, v = frames[i]
        ring_offset.append(len(pts))
        
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
        for k in range(n_circ):
            a = r0 + k
            b = r0 + ((k + 1) % n_circ)
            c_idx = r1 + ((k + 1) % n_circ)
            d_idx = r1 + k
            quads.append((a, b, c_idx, d_idx))
    
    return pts, quads, modal_radii
