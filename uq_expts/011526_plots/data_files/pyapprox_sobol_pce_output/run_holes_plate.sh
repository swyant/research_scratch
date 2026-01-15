#!/bin/bash
# run_holes_plate.sh
# Reads params.in, substitutes into template, runs holes_plate, extracts results

# Paths (must be absolute or relative from work directory)
HOLES_PLATE_EXE="/home/swyant/local_software/sumMITApps/Dakota_apps/holes_in_plate/build/holes_plate"
MESH_FILE="hole.summit"
TEMPLATE_FILE="holeplate.template"
#PYTHON_EXE="/home/swyant/cesmix/exploratory/private_research/uq_expts/summit_tests/pyapprox_expts/.venv/bin"
PYTHON_EXE="python"

# Read parameters from params.in (one value per line: rho, E, nu)
RHO=$(sed -n '1p' params.in)
E=$(sed -n '2p' params.in)
NU=$(sed -n '3p' params.in)

# Substitute parameters into template to create input.dat
# Template format: "1\n13 elastic {rho=1.0} {E=1.0} {nu=0.27} 0\n"
# Note: holes_plate expects input.dat (not input.in) in the current directory
sed -e "s/{rho=1.0}/${RHO}/" \
    -e "s/{E=1.0}/${E}/" \
    -e "s/{nu=0.27}/${NU}/" \
    "${TEMPLATE_FILE}" > material.dat
#    "${TEMPLATE_FILE}" > input.dat

# Run holes_plate with the mesh file
# The executable expects:
#   - Mesh file as argument
#   - input.dat in current directory
#   - Outputs to results.csv
"${HOLES_PLATE_EXE}" "${MESH_FILE}" > simulation_output.txt 2>&1

# Check if execution was successful
if [ ! -f "results.csv" ]; then
    echo "ERROR: holes_plate did not generate results.csv" > results.out
    echo "Stdout/Stderr:" >> results.out
    cat simulation_output.txt >> results.out
    echo "Params: rho=$RHO, E=$E, nu=$NU" >> results.out
    exit 1
fi

# Parse results.csv and write to results.out (4 QoI values, one per line)
# CSV format: E_effective,stress,strain,max_von_Mises_stress
"${PYTHON_EXE}" << PYTHON_EOF
import csv
import numpy as np

results = []

try:
    with open('results.csv', 'r') as f:
        reader = csv.DictReader(f)
        row = next(reader)
        
        # Extract QoI in the correct order
        results.append(float(row['E_effective']))
        results.append(float(row['stress']))
        results.append(float(row['max_von_Mises_stress']))
    
    # Write one value per line (AynchModel expects this format with default np.loadtxt)
    np.savetxt('results.out', results)
    
except Exception as e:
    print(f"Error parsing results.csv: {e}")
    np.savetxt('results.out', [np.nan, np.nan, np.nan])
    exit(1)
PYTHON_EOF
