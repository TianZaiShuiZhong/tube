from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

from ..model import Model, FourierModalState


def _safe_dir(p0: tuple[float, float, float], p1: tuple[float, float, float]):
    dx = (p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2])
    n = math.sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2])
    if n == 0:
        return (1.0, 0.0, 0.0)
    return (dx[0] / n, dx[1] / n, dx[2] / n)


def _orthonormal(d: tuple[float, float, float]):
    # build two orthonormal vectors perpendicular to d
    # choose an arbitrary vector not parallel to d
    ax = (0.0, 0.0, 1.0) if abs(d[2]) < 0.9 else (0.0, 1.0, 0.0)
    # u = d x ax
    u = (
        d[1] * ax[2] - d[2] * ax[1],
        d[2] * ax[0] - d[0] * ax[2],
        d[0] * ax[1] - d[1] * ax[0],
    )
    nu = math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
    if nu == 0:
        u = (1.0, 0.0, 0.0)
    else:
        u = (u[0] / nu, u[1] / nu, u[2] / nu)
    # v = d x u
    v = (
        d[1] * u[2] - d[2] * u[1],
        d[2] * u[0] - d[0] * u[2],
        d[0] * u[1] - d[1] * u[0],
    )
    return u, v


    return u, v


def _extract_paths(model: Model) -> list[list[int]]:
    """Decompose the model into a list of connected node paths (polylines)."""
    adj: dict[int, list[int]] = {}
    for e in model.elements.values():
        n1, n2 = e.end_nodes()
        adj.setdefault(n1, []).append(n2)
        adj.setdefault(n2, []).append(n1)

    # Find start nodes (degree != 2)
    start_nodes = [n for n, neighbors in adj.items() if len(neighbors) != 2]
    # If no start nodes but elements exist, it's a ring
    if not start_nodes and adj:
        start_nodes = [next(iter(adj))]

    visited_edges: set[tuple[int, int]] = set()
    paths: list[list[int]] = []

    def get_edge_key(u, v):
        return tuple(sorted((u, v)))

    for start_node in start_nodes:
        # Explore all branches from this node
        neighbors = adj.get(start_node, [])
        for next_node in neighbors:
            edge = get_edge_key(start_node, next_node)
            if edge in visited_edges:
                continue
            
            # Start a new path
            current_path = [start_node, next_node]
            visited_edges.add(edge)
            
            curr = next_node
            while True:
                # Find next step
                # Look for neighbors of curr excluding the one we came from (current_path[-2])
                prev = current_path[-2]
                candidates = [n for n in adj[curr] if n != prev]
                
                if len(candidates) != 1:
                    # Junction or endpoint -> stop path
                    break
                
                next_step = candidates[0]
                edge_next = get_edge_key(curr, next_step)
                if edge_next in visited_edges:
                    break # Should not happen in tree, but safety for loops
                
                visited_edges.add(edge_next)
                current_path.append(next_step)
                curr = next_step
            
            paths.append(current_path)
            
    # Check for any unvisited edges (disconnected rings)
    # Naive check: iterate all edges, if not visited, start a path
    for e in model.elements.values():
        n1, n2 = e.end_nodes()
        if get_edge_key(n1, n2) not in visited_edges:
            # Found a disconnected component (likely a ring)
            # Trace it
            curr = n1
            path = [curr]
            while True:
                # Find unvisited neighbor
                found_next = False
                for nxt in adj[curr]:
                    ek = get_edge_key(curr, nxt)
                    if ek not in visited_edges:
                        visited_edges.add(ek)
                        path.append(nxt)
                        curr = nxt
                        found_next = True
                        break
                if not found_next:
                    break
            paths.append(path)

    return paths


def write_pipe_surface_quads(
    path: Path,
    model: Model,
    centerline_nodes: Sequence[int] | None,
    radius: float,
    n_circ: int = 24,
    modal_state: FourierModalState | None = None,
):
    """Write pipe surface with optional modal displacement field for visualization."""
    # Import here to avoid circular import
    from ..solver.modal import expand_surface_with_modal
    
    # Decompose model into paths to handle branching and ordering correctly
    paths = _extract_paths(model)
    
    all_pts = []
    all_quads = []
    all_modal_radii = []
    
    for p_nodes in paths:
        if len(p_nodes) < 2:
            continue
            
        # Generate surface for this path segment
        if modal_state is not None and modal_state.modal_amplitudes:
            pts, quads, modal_radii = expand_surface_with_modal(
                model, p_nodes, radius, n_circ, modal_state
            )
        else:
            # Original simple expansion (no modal)
            pts = []
            quads = []
            modal_radii = []
            
            # Build center coordinates, removing consecutive duplicates (zero-length segments)
            raw_centers = [model.nodes[nid] for nid in p_nodes]
            centers: list = []
            node_map: list = []
            tol = 1e-9
            for nd in raw_centers:
                if not centers:
                    centers.append(nd)
                    node_map.append(nd)
                else:
                    prev = centers[-1]
                    if abs(prev.x - nd.x) < tol and abs(prev.y - nd.y) < tol and abs(prev.z - nd.z) < tol:
                        # skip duplicate center
                        continue
                    centers.append(nd)
                    node_map.append(nd)

            # Compute smooth tangents using forward/backward averaging to avoid flipped frames
            frames = []
            tangents = []
            for i in range(len(centers)):
                if i == 0:
                    p0 = (centers[i].x, centers[i].y, centers[i].z)
                    p1 = (centers[i + 1].x, centers[i + 1].y, centers[i + 1].z)
                    d1 = _safe_dir(p0, p1)
                    d = d1
                elif i == len(centers) - 1:
                    p0 = (centers[i - 1].x, centers[i - 1].y, centers[i - 1].z)
                    p1 = (centers[i].x, centers[i].y, centers[i].z)
                    d2 = _safe_dir(p0, p1)
                    d = d2
                else:
                    p_prev = (centers[i - 1].x, centers[i - 1].y, centers[i - 1].z)
                    p_curr = (centers[i].x, centers[i].y, centers[i].z)
                    p_next = (centers[i + 1].x, centers[i + 1].y, centers[i + 1].z)
                    v_prev = (p_curr[0] - p_prev[0], p_curr[1] - p_prev[1], p_curr[2] - p_prev[2])
                    v_next = (p_next[0] - p_curr[0], p_next[1] - p_curr[1], p_next[2] - p_curr[2])
                    d_prev = _safe_dir(p_prev, p_curr)
                    d_next = _safe_dir(p_curr, p_next)
                    # Use a length-weighted bisector to round corners smoothly
                    len_prev = math.sqrt(v_prev[0] * v_prev[0] + v_prev[1] * v_prev[1] + v_prev[2] * v_prev[2])
                    len_next = math.sqrt(v_next[0] * v_next[0] + v_next[1] * v_next[1] + v_next[2] * v_next[2])
                    w_prev = len_prev if len_prev > 0 else 1.0
                    w_next = len_next if len_next > 0 else 1.0
                    d = (
                        w_prev * d_prev[0] + w_next * d_next[0],
                        w_prev * d_prev[1] + w_next * d_next[1],
                        w_prev * d_prev[2] + w_next * d_next[2],
                    )
                    nd = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                    if nd == 0:
                        d = d_next
                    else:
                        d = (d[0] / nd, d[1] / nd, d[2] / nd)
                tangents.append(d)
                u, v = _orthonormal(d)
                # ensure continuity: flip if inconsistent with previous u
                if frames:
                    pu, pv = frames[-1]
                    dot = pu[0] * u[0] + pu[1] * u[1] + pu[2] * u[2]
                    if dot < 0:
                        u = (-u[0], -u[1], -u[2])
                        v = (-v[0], -v[1], -v[2])
                frames.append((u, v))
            
            ring_offset = []
            for i, c in enumerate(centers):
                u, v = frames[i]
                ring_offset.append(len(pts))
                # per-node radius: if different radii are expected per section, user can extend here
                r = radius
                for k in range(n_circ):
                    theta = 2.0 * math.pi * (k / n_circ)
                    px = c.x + r * (math.cos(theta) * u[0] + math.sin(theta) * v[0])
                    py = c.y + r * (math.cos(theta) * u[1] + math.sin(theta) * v[1])
                    pz = c.z + r * (math.cos(theta) * u[2] + math.sin(theta) * v[2])
                    pts.append((px, py, pz))
                    modal_radii.append(0.0)
            
            for i in range(len(centers) - 1):
                r0 = ring_offset[i]
                r1 = ring_offset[i + 1]
                for k in range(n_circ):
                    a = r0 + k
                    b = r0 + ((k + 1) % n_circ)
                    c_idx = r1 + ((k + 1) % n_circ)
                    d_idx = r1 + k
                    quads.append((a, b, c_idx, d_idx))

        # Merge into global list
        base_idx = len(all_pts)
        all_pts.extend(pts)
        all_modal_radii.extend(modal_radii)
        for (a, b, c, d) in quads:
            all_quads.append((a + base_idx, b + base_idx, c + base_idx, d + base_idx))
    
    # Write legacy VTK POLYDATA with quads and modal scalars
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("tubo pipe surface with modal\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(all_pts)} float\n")
        for x, y, z in all_pts:
            f.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        nquads = len(all_quads)
        f.write(f"POLYGONS {nquads} {nquads * 5}\n")
        for a, b, c, d in all_quads:
            f.write(f"4 {a} {b} {c} {d}\n")
        
        # Write modal radial component as scalar for visualization
        if all_modal_radii:
            f.write(f"POINT_DATA {len(all_pts)}\n")
            f.write("SCALARS modal_radial_disp float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for val in all_modal_radii:
                f.write(f"{val:.9g}\n")
