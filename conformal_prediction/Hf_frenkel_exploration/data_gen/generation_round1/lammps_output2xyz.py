#!/usr/bin/env python 
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

# only needed because lammpsrun.py in ASE has the timestep commmented out right now
def extract_timesteps(dump_file):
    """Extract timesteps from the LAMMPS dump file."""
    timesteps = []
    with open(dump_file, 'r') as file:
        for line in file:
            if line.startswith('ITEM: TIMESTEP'):
                timestep = int(next(file).strip())
                timesteps.append(timestep)
    return timesteps

def parse_lammpsout_pe(out_file):
    """
    Extracts the Step and PotEng columns from a LAMMPS output file
    
    Args:
        output_file (str): Path to the LAMMPS output file
        
    Returns:
        dict: Dictionary with Step values as keys and PotEng values as values
    """
    result = {}
    data_section = False
    column_indices = {}
    
    with open(out_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Check if we're at the header line
            if "Step" in line and "PotEng" in line:
                # Parse the header to find column indices
                headers = line.split()
                for i, header in enumerate(headers):
                    column_indices[header] = i
                data_section = True
                continue
                
            # If we're in the data section and line starts with a number
            if data_section and line[0].isdigit():
                # Check if we've reached a line that doesn't start with numbers
                if not line.split()[0].isdigit():
                    data_section = False
                    continue
                    
                parts = line.split()
                # Make sure we have enough columns
                if len(parts) > column_indices["PotEng"]:
                    step = int(parts[column_indices["Step"]])
                    poteng = float(parts[column_indices["PotEng"]])
                    result[step] = poteng
            
            # If we encounter a line starting with "Loop time", we're done with data
            if data_section and line.startswith("Loop time"):
                break
                
    return result


def main(out_file, dump_file, output_extxyz):
    # Parse the output xyz file
    timestep_energy = parse_lammpsout_pe(out_file)

    # Read configurations from the LAMMPS dump file
    configurations = read(dump_file, format='lammps-dump-text', index=':')

    # Extract timesteps from the dump file
    timesteps = extract_timesteps(dump_file)

    # Ensure the number of timesteps matches the number of configurations
    if len(timesteps) != len(configurations):
        raise ValueError("Mismatch between number of timesteps and configurations.")

    # List to hold Atoms objects with energies
    atoms_list = []

    for atoms,timestep in zip(configurations,timesteps):
        #timestep = atoms.info.get('timestep')
        if timestep in timestep_energy:
            # Assign the potential energy to the Atoms object
            #atoms.info['energy'] = timestep_energy[timestep]
            #atoms_list.append(atoms)
            existing_forces = atoms.get_forces()
            atoms.calc = SinglePointCalculator(atoms, energy=timestep_energy[timestep],forces=existing_forces)
            atoms_list.append(atoms)
        else:
            print(f"Warning: No energy found for timestep {timestep}. Skipping this configuration.")

    # Write the configurations to an extended XYZ file
    if atoms_list:
        write(output_extxyz, atoms_list, format='extxyz')
        print(f"Successfully wrote {len(atoms_list)} configurations to {output_extxyz}.")
        
        # Write the configurations to a plain XYZ file without energies and forces
        # GPT changed the format here to just 'xyz', but you want to keep it 'extxyz' to pass plain 
        # actually don't want to pass plain=True, because lose the lattice 
        #write(output_xyz, atoms_list, format='extxyz', write_results=False)
        #print(f"Successfully wrote {len(atoms_list)} configurations to {output_xyz}.")
    else:
        print("No configurations were written. Please check the input files.")

if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) != 4:
        print("huh")
        sys.exit(1)

    outdir = sys.argv[1]
    dump_fname = sys.argv[2]
    extxyz_fname = sys.argv[3]

    lammps_outfile = os.path.join(outdir, 'out.txt')
    dump_file = os.path.join(outdir,dump_fname)
   
    main(lammps_outfile,dump_file,extxyz_fname)
