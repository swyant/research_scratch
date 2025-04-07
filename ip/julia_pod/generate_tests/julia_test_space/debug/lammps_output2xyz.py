#!/usr/bin/env python 
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

def parse_pe_file(pe_file):
    """Parse the pe.dat file to extract timestep-energy pairs."""
    timestep_energy = {}
    with open(pe_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.split()
            timestep = int(parts[0])
            energy = float(parts[1])
            timestep_energy[timestep] = energy
    return timestep_energy

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

def main(dump_file, pe_file, output_extxyz, output_xyz):
    # Parse the potential energy file
    timestep_energy = parse_pe_file(pe_file)

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
        write(output_xyz, atoms_list, format='extxyz', write_results=False)
        print(f"Successfully wrote {len(atoms_list)} configurations to {output_xyz}.")
    else:
        print("No configurations were written. Please check the input files.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print("Usage: python script.py <dump_file> <pe_file> <output_file>")
        sys.exit(1)
    dump_file = sys.argv[1]
    pe_file = sys.argv[2]
    output_extxyz = sys.argv[3]
    output_xyz = sys.argv[4]
    main(dump_file, pe_file, output_extxyz, output_xyz)

'''
chatgpt free 3/27/25

I have two lammps output files, which I want to take and convert to an extxyz file 
The first is a dump.custom, a snippet of this file is provide here:
"""
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
ITEM: ATOMS id type x y z fx fy fz
   1 1    0.00000000000000000    0.00000000000000000    0.00000000000000000    0.00000005120123792    0.00000005120122788    0.00000005120125162
   2 1    1.51609879000000003    1.51609879000000003    1.51609879000000003   -0.00000005120123819   -0.00000005120122869   -0.00000005120125253
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
ITEM: ATOMS id type x y z fx fy fz
   1 1    3.00593846161599210    3.00593846161599254    3.00593846161599254    1.06475587551207296    1.06475587551207540    1.06475587551207940
   2 1    1.54235790332559097    1.54235790332559097    1.54235790332559053   -1.06475587551207296   -1.06475587551207562   -1.06475587551207962
ITEM: TIMESTEP
200
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
0.0000000000000000e+00 3.0321975749415837e+00
ITEM: ATOMS id type x y z fx fy fz
   1 1    3.03211250138699251    3.03211250138699251    3.03211250138699340    0.00344450244792491    0.00344450244792023    0.00344450244787013
   2 1    1.51618386355459100    1.51618386355459100    1.51618386355458989   -0.00344450244792460   -0.00344450244791932   -0.00344450244787015
"""
A second file is pe.dat, which outputs the potential energy every 100 timesteps:
"""
# Fix print output for fix eprint
0 -10.497338278136052026
100 -10.413519428852257676
200 -10.497337399014813997
300 -10.241507666893017614
400 -10.497168958758948065
500 -10.381284690476025645
"""
Reading in these two files, I want to output an extxyz file with the corresponding snapshots, including both the forces (obtained from the dump.custom) and the potential energies (from pe.dat). Please leverage the atomic simulation environment, and please ensure that the timesteps line up between snapshots in the dump.custom and the entries in pe.dat

...

How can I modify the above script to write a second xyz file that doesn't contain any information about either forces or energies
'''
