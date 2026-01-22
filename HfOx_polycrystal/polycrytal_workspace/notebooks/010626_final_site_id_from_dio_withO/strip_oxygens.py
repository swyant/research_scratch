#!/usr/bin/env python3
"""
Script to remove oxygen atoms from a structure file and write to LAMMPS DATA format.
Reads xyz or lammps-dump-text files and optionally excludes a single oxygen atom.
"""

import argparse
from ase.io import read, write

from pathlib import Path
import json

def main():
    parser = argparse.ArgumentParser(
        description='Remove oxygen atoms from structure file and write to LAMMPS DATA format'
    )
    parser.add_argument(
        '--exclude',
        type=int,
        default=None,
        help='Atom ID (0-indexed) to exclude from removal'
    )
    parser.add_argument(
        'input_file',
        help='Input file (xyz or lammps-dump-text)'
    )
    parser.add_argument(
        'output_file',
        help='Output LAMMPS DATA file'
    )
    
    args = parser.parse_args()
    
    # Read the structure file
    atoms = read(args.input_file)
    
    # Create a list of indices to keep
    indices_to_keep = []
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'O':
            # Keep this oxygen only if it matches the excluded ID
            if args.exclude is not None and i == args.exclude:
                indices_to_keep.append(i)
        else:
            # Keep all non-oxygen atoms
            indices_to_keep.append(i)
    
    # Create new atoms object with only the atoms we want to keep
    filtered_atoms = atoms[indices_to_keep]

    ########### DIAGNOSTIC ##############
    hf_index_map = {}
    new_index = 0
    
    for original_index in indices_to_keep:
        if atoms[original_index].symbol == 'Hf':
            hf_index_map[new_index] = original_index
        new_index += 1     

    mapping_file = Path(args.output_file).with_suffix('.hf_mapping.json')
    with open(mapping_file, 'w') as f:
        json.dump(hf_index_map, f, indent=2)
    ########### DIAGNOSTIC ##############

    # Write to LAMMPS DATA format
    write(
        args.output_file,
        filtered_atoms,
        format='lammps-data',
        specorder=['Hf', 'O'],
        masses=True
    )
    
    print(f"Removed {len(atoms) - len(filtered_atoms)} oxygen atom(s)")
    print(f"Kept {len(filtered_atoms)} atom(s)")
    print(f"Written to {args.output_file}")

if __name__ == '__main__':
    main()
