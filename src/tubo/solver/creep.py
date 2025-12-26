from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..model import Model
from .assembly import assemble_linear_system, ElementKinematics


@dataclass
class IntegrationPoint:
    y: float
    z: float
    weight: float
    eps_cr: float = 0.0  # accumulated creep strain
    eps_pl: float = 0.0  # accumulated plastic strain
    eps_pl_eq: float = 0.0  # equivalent plastic strain (for hardening)


def _generate_integration_points(ro: float, ri: float, n_circ: int = 12, n_rad: int = 2) -> list[IntegrationPoint]:
    points = []
    dr = (ro - ri) / n_rad
    d_theta = 2.0 * math.pi / n_circ
    
    for i in range(n_rad):
        r_mid = ri + (i + 0.5) * dr
        # Area of the ring sector
        area_ring = math.pi * ((r_mid + 0.5*dr)**2 - (r_mid - 0.5*dr)**2)
        weight = area_ring / n_circ
        
        for j in range(n_circ):
            theta = j * d_theta
            y = r_mid * math.cos(theta)
            z = r_mid * math.sin(theta)
            points.append(IntegrationPoint(y, z, weight))
            
    return points


def solve_creep_time_series(model: Model, *, time_end: float, nsubsteps: int, enable_pressure_equiv: bool = True):
    # High-efficiency creep solver with cross-section integration
    if not model.materials:
        raise ValueError("No materials")
    mat = model.materials[min(model.materials.keys())]
    
    # Creep parameters (optional)
    c1 = c2 = c3 = c4 = 0.0
    if mat.creep:
        c1, c2, c3, c4 = mat.creep
    
    # Plasticity parameters
    sy, etan = 1e20, 0.0  # Default infinite yield
    if mat.plasticity:
        sy, etan = mat.plasticity
    
    # Section properties
    if not model.sections:
        raise ValueError("No sections")
    sec = model.sections[min(model.sections.keys())]
    ro = 0.5 * sec.outer_diameter
    ri = ro - sec.thickness
    
    # Initialize integration points for each element
    # Use 12 circ, 2 rad (High-efficiency: enough for bending and thickness effect)
    ips_template = _generate_integration_points(ro, ri, n_circ=12, n_rad=2)
    
    # Deep copy for each element
    elem_ips: dict[int, list[IntegrationPoint]] = {}
    for eid in model.elements:
        elem_ips[eid] = [IntegrationPoint(ip.y, ip.z, ip.weight) for ip in ips_template]
        
    times = np.linspace(0.0, float(time_end), int(nsubsteps))
    dt = np.diff(times, prepend=0.0)
    
    series = []
    
    # Initial state (creep strains)
    eps0_by_elem: dict[int, float] = {}
    curv_by_elem: dict[int, tuple[float, float]] = {}
    
    for t_idx, (t, dti) in enumerate(zip(times, dt)):
        # 1. Solve elastic with current initial strains
        node_ids, node_index, K, F, elem_infos, dof_map, props = assemble_linear_system(
            model, 
            extra_initial_strain_by_elem=eps0_by_elem,
            extra_initial_curvature_by_elem=curv_by_elem,
            enable_pressure_equiv=enable_pressure_equiv
        )
        E, A, _ = props
        
        # Solve system
        ndof = K.shape[0]
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
                modal_amplitudes[nid] = {2: float(amp)}
            
            from ..model import FourierModalState
            model.modal_state = FourierModalState(modal_amplitudes=modal_amplitudes)

        series.append((float(t), {"node_ids": np.array(node_ids, dtype=int), "U": U}))
        
        # 2. Update creep and plastic strains
        # Note: For t=0 (dti=0), we still check plasticity (instantaneous)
            
        for info in elem_infos:
            # Recover element displacements
            u_global = U[info.dofs]
            u_local = info.T @ u_global
            
            # Calculate generalized strains at midpoint
            eps_ax = (u_local[6] - u_local[0]) / info.L
            
            # Curvatures at midpoint
            # kz (bending in xy plane)
            kz_mid = (u_local[11] - u_local[5]) / info.L
            # ky (bending in xz plane)
            ky_mid = (u_local[10] - u_local[4]) / info.L
            
            # Update integration points
            ips = elem_ips[info.eid]
            
            sum_eps_inelastic = 0.0
            sum_eps_inelastic_y = 0.0 
            sum_eps_inelastic_z = 0.0 
            
            total_area = 0.0
            total_Iy = 0.0
            total_Iz = 0.0
            
            for ip in ips:
                # Total strain at point
                # eps = eps0 - y*kz + z*ky
                eps_total = eps_ax - ip.y * kz_mid + ip.z * ky_mid
                
                # Trial Stress (elastic prediction)
                # sigma_trial = E * (eps_total - eps_cr - eps_pl_old)
                sigma_trial = E * (eps_total - ip.eps_cr - ip.eps_pl)
                
                # --- Plasticity Return Mapping (1D) ---
                # Yield function: f = |sigma| - (sy + H * eps_pl_eq)
                # Hardening modulus H = E_tan * E / (E - E_tan)
                H = 0.0
                if etan < E:
                    H = etan * E / (E - etan)
                
                current_yield = sy + H * ip.eps_pl_eq
                f_yield = abs(sigma_trial) - current_yield
                
                sigma_final = sigma_trial
                
                if f_yield > 0:
                    # Plastic correction
                    # d_lambda = f / (E + H)
                    d_lambda = f_yield / (E + H)
                    sign_sigma = 1.0 if sigma_trial >= 0 else -1.0
                    
                    # Update plastic strain
                    d_eps_pl = d_lambda * sign_sigma
                    ip.eps_pl += d_eps_pl
                    ip.eps_pl_eq += d_lambda
                    
                    # Return to yield surface
                    sigma_final = sigma_trial - E * d_eps_pl
                
                # --- Creep Update ---
                if dti > 0:
                    # Creep rate based on final stress state
                    # Norton: B * sigma^n
                    s_val = abs(sigma_final)
                    s_sign = 1.0 if sigma_final >= 0 else -1.0
                    
                    rate = c1 * (s_val ** c2) * s_sign
                    d_eps_cr = rate * dti
                    ip.eps_cr += d_eps_cr
                
                # Accumulate total inelastic strain for next step's load vector
                eps_inelastic = ip.eps_cr + ip.eps_pl
                
                # Integrate
                sum_eps_inelastic += eps_inelastic * ip.weight
                sum_eps_inelastic_y += eps_inelastic * ip.y * ip.weight
                sum_eps_inelastic_z += eps_inelastic * ip.z * ip.weight
                
                total_area += ip.weight
                total_Iy += (ip.z ** 2) * ip.weight
                total_Iz += (ip.y ** 2) * ip.weight
            
            # Equivalent initial generalized strains
            if total_area > 0:
                eps0_by_elem[info.eid] = sum_eps_inelastic / total_area
            
            if total_Iz > 0:
                kz_cr = -sum_eps_inelastic_y / total_Iz
            else:
                kz_cr = 0.0
                
            if total_Iy > 0:
                ky_cr = sum_eps_inelastic_z / total_Iy
            else:
                ky_cr = 0.0
            
            curv_by_elem[info.eid] = (ky_cr, kz_cr)
            
    return series
