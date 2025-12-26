
import sys
import math
import re

def read_vtk_max_disp(filepath):
    """Reads a legacy VTK file and returns the maximum displacement magnitude."""
    max_disp = 0.0
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        # Find VECTORS displacement
        start_idx = -1
        for i, line in enumerate(lines):
            if line.startswith("VECTORS displacement"):
                start_idx = i + 1
                break
        
        if start_idx == -1:
            print(f"Error: Could not find displacement vectors in {filepath}")
            return 0.0
            
        # Read vectors until we hit another section or EOF
        for i in range(start_idx, len(lines)):
            line = lines[i].strip()
            if not line or not line[0].replace('-','').replace('.','').isdigit():
                # Stop if we hit a non-number line (like SCALARS)
                # But wait, SCALARS line starts with SCALARS
                if line.startswith("SCALARS") or line.startswith("LOOKUP_TABLE"):
                    break
            
            try:
                parts = line.split()
                if len(parts) == 3:
                    ux, uy, uz = map(float, parts)
                    mag = math.sqrt(ux*ux + uy*uy + uz*uz)
                    if mag > max_disp:
                        max_disp = mag
            except ValueError:
                break
                
    except FileNotFoundError:
        print(f"Error: File not found {filepath}")
        return 0.0
        
    return max_disp

def verify_pipe288_plast():
    print("\n--- Verifying PIPE288_PLAST ---")
    disp = read_vtk_max_disp("build/out_pipe288_plast.vtk")
    print(f"Max Displacement (FEM): {disp:.4f} mm")
    
    # Theoretical Estimate
    # Thermal (200C, alpha=1.2e-5, L=2500) -> ~3.15 mm (Vector sum of segments)
    # Gravity (Cantilever approx) -> ~0.7 mm (Axial + Bending)
    # Total vector sum approx 3.5 mm
    
    expected = 3.5
    print(f"Theoretical Estimate: ~{expected:.2f} mm")
    
    if abs(disp - expected) < 1.0: # Allow some margin for BCs and exact gravity integration
        print("✅ Result is within reasonable range of theory.")
    else:
        print("⚠️ Result deviates significantly from rough theory. Check loads.")

def verify_elbow290_plast():
    print("\n--- Verifying ELBOW290_PLAST ---")
    disp = read_vtk_max_disp("build/out_elbow290_plast.vtk")
    print(f"Max Displacement (FEM): {disp:.4f} mm")
    
    # Thermal (350C, alpha=1.2e-5, L~1047) -> ~4.08 mm
    # Gravity -> ?
    # Moment -> ?
    # Just checking if it ran and produced deformation
    if disp > 0.1:
        print("✅ Significant deformation detected.")
    else:
        print("⚠️ Deformation is suspiciously small.")

def verify_pipe288_creep():
    print("\n--- Verifying PIPE288_CREEP ---")
    d0 = read_vtk_max_disp("build/pipe288_creep/series_0000.vtk")
    d_end = read_vtk_max_disp("build/pipe288_creep/series_0099.vtk") # Assuming 100 steps
    
    print(f"Initial Displacement (t=0): {d0:.4f} mm")
    print(f"Final Displacement (t=end): {d_end:.4f} mm")
    
    if abs(d_end - d0) > 0.1:
        print("✅ Creep deformation detected (Displacement changed over time).")
    else:
        print("⚠️ No significant creep deformation detected.")

if __name__ == "__main__":
    verify_pipe288_plast()
    verify_elbow290_plast()
    verify_pipe288_creep()
