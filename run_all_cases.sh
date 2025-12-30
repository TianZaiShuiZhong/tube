#!/bin/bash
set -e

# Cleanup
echo "Cleaning up old outputs..."
rm -rf build out_*.vtk pipe288_creep* elbow290_plast*
mkdir -p build

# 1. PIPE288_PLAST
echo "Running PIPE288_PLAST..."
PYTHONPATH=src python3 -m tubo run --cdb "开放原子-管单元/PIPE288_PLAST.cdb" --out "build/out_pipe288_plast.vtk" --surface "build/out_pipe288_plast_surface.vtk" --pressure-equiv

# 2. ELBOW290_PLAST
echo "Running ELBOW290_PLAST..."
PYTHONPATH=src python3 -m tubo run --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" --out "build/out_elbow290_plast.vtk" --surface "build/out_elbow290_plast_surface.vtk" --pressure-equiv

# 3. PIPEMIX_PLAST
echo "Running PIPEMIX_PLAST..."
PYTHONPATH=src python3 -m tubo run --cdb "开放原子-管单元/PIPEMIX_PLAST.cdb" --out "build/out_pipemix_plast.vtk" --surface "build/out_pipemix_plast_surface.vtk" --pressure-equiv

# 4. PIPE288_CREEP
echo "Running PIPE288_CREEP..."
PYTHONPATH=src python3 -m tubo creep --cdb "开放原子-管单元/PIPE288_CREEP.cdb" --outdir "build/pipe288_creep" --basename "series" --pressure-equiv

# 5. ELBOW290_CREEP
echo "Running ELBOW290_CREEP..."
PYTHONPATH=src python3 -m tubo creep --cdb "开放原子-管单元/ELBOW290_CREEP.cdb" --outdir "build/elbow290_creep" --basename "series" --pressure-equiv

# 6. PIPEMIX_CREEP
echo "Running PIPEMIX_CREEP..."
PYTHONPATH=src python3 -m tubo creep --cdb "开放原子-管单元/PIPEMIX_CREEP.cdb" --outdir "build/pipemix_creep" --basename "series" --pressure-equiv

echo "All cases completed."
