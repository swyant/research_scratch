#!/usr/bin/env python3
"""
Script to convert VASP vasprun.xml files to concatenated XYZ files.
Creates two versions: one with standard ASE fields and one with REF_ prefixed fields.
"""

import sys
from pathlib import Path
from ase.io import read, write
from ase import Atoms
import numpy as np
import copy 

def voigt_6_to_full_3x3_stress(stress_vector):
    """
    Form a 3x3 stress matrix from a 6 component vector in Voigt notation
    """
    s1, s2, s3, s4, s5, s6 = stress_vector
    return np.array([[s1, s6, s5],
                     [s6, s2, s4],
                     [s5, s4, s3]])


def parse_vasprun_to_atoms(vasprun_path):
    """
    Parse a vasprun.xml file and return the final ASE Atoms object.
    
    Args:
        vasprun_path: Path to vasprun.xml file
        
    Returns:
        Final Atoms object from the calculation, or None if error
    """
    try:
        # Read only the last structure (index=-1)
        atoms = read(vasprun_path, index=-1, format='vasp-xml')
        return atoms
    except Exception as e:
        print(f"Error reading {vasprun_path}: {e}", file=sys.stderr)
        return None


def create_standard_atoms(atoms):
    """
    Create atoms object with standard ASE fields.
    Forces, stress, energy, and free_energy in standard locations.
    
    Args:
        atoms: ASE Atoms object from VASP
        
    Returns:
        ASE Atoms object with standard fields
    """
    new_atoms = copy.deepcopy(atoms)
    
    # Get calculator results if available
    if atoms.calc is not None and hasattr(atoms.calc, 'results'):
        results = atoms.calc.results
        
        # Set forces (standard ASE array)
        if 'forces' in results:
            new_atoms.arrays['forces'] = results['forces'].copy()
        
        # Create a simple calculator to store results
        from ase.calculators.singlepoint import SinglePointCalculator
        
        calc_results = {}
       
         # Store energy and free_energy (should be the same)
        if 'free_energy' in results:
            energy = results['free_energy']
            calc_results['energy'] = energy
            calc_results['free_energy'] = energy


        # Store stress in standard location
        if 'stress' in results:
            calc_results['stress'] = results['stress'].copy()
         
        # Attach calculator with results
        calc = SinglePointCalculator(new_atoms, **calc_results)
        new_atoms.calc = calc
    
    return new_atoms


def create_ref_atoms(atoms):
    """
    Create atoms object with REF_ prefixed fields.
    Forces -> REF_forces, stress -> REF_stress, energy -> REF_energy
    
    Args:
        atoms: ASE Atoms object from VASP
        
    Returns:
        ASE Atoms object with REF_ prefixed fields
    """
    new_atoms = copy.deepcopy(atoms)
    
    # Get calculator results if available
    if atoms.calc is not None and hasattr(atoms.calc, 'results'):
        results = atoms.calc.results
        
        # Store forces in REF_forces array
        if 'forces' in results:
            new_atoms.arrays['REF_forces'] = results['forces'].copy()
            # Remove standard forces array if it exists
            if 'forces' in new_atoms.arrays:
                del new_atoms.arrays['forces']
        
        # Create calculator to store REF_ results
        from ase.calculators.singlepoint import SinglePointCalculator
        
        calc_results = {}
        atoms_info = {} 
        
         
        # Store stress as REF_stress
        if 'stress' in results:
            stress_voigt = results['stress']
            stress_3x3 = voigt_6_to_full_3x3_stress(stress_voigt)           
            atoms_info['REF_stress'] = stress_3x3.flatten()
        
        # Store energy as REF_energy and keep free_energy
        if 'free_energy' in results:
            energy = results['free_energy']
            calc_results['free_energy'] = energy
            atoms_info['REF_energy'] = energy
        
        # Attach calculator with results
        calc = SinglePointCalculator(new_atoms, **calc_results)
        new_atoms.calc = calc
        new_atoms.info = atoms_info
    
    return new_atoms


def process_file_list(input_file, output_standard, output_ref):
    """
    Process list of vasprun.xml files and create two concatenated XYZ outputs.
    
    Args:
        input_file: Path to text file containing list of vasprun.xml paths
        output_standard: Path for output XYZ with standard fields
        output_ref: Path for output XYZ with REF_ fields
    """
    # Read the list of file paths
    with open(input_file, 'r') as f:
        file_paths = [line.strip() for line in f if line.strip()]
    
    print(f"Processing {len(file_paths)} files...")
    
    standard_atoms_list = []
    ref_atoms_list = []
    
    processed = 0
    skipped = 0
    
    for path_str in file_paths:
        path = Path(path_str)
        
        if not path.exists():
            print(f"Warning: File not found: {path}", file=sys.stderr)
            skipped += 1
            continue
        
        # Parse vasprun.xml
        atoms = parse_vasprun_to_atoms(path)
        
        if not atoms:
            print(f"{path} skipped")
            skipped += 1
            continue
        
        # Create both versions
        standard_atoms = create_standard_atoms(atoms)
        #standard_atoms = copy.deepcopy(atoms)
        ref_atoms = create_ref_atoms(atoms)
        
        standard_atoms_list.append(standard_atoms)
        ref_atoms_list.append(ref_atoms)
        
        processed += 1
        if processed % 10 == 0:
            print(f"Processed {processed}/{len(file_paths)} files...")
    
    print(f"\nTotal structures: {len(standard_atoms_list)}")
    print(f"Files processed: {processed}")
    print(f"Files skipped: {skipped}")
    
    # Write concatenated XYZ files
    if standard_atoms_list:
        print(f"\nWriting standard format to: {output_standard}")
        write(output_standard, standard_atoms_list, format='extxyz')
        
        print(f"Writing REF format to: {output_ref}")
        write(output_ref, ref_atoms_list, format='extxyz')
        
        print("\nDone!")
    else:
        print("\nNo structures to write!", file=sys.stderr)
        sys.exit(1)


def main():
    """Main entry point."""
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file_list> <output_standard.xyz> <output_ref.xyz>")
        print("\nExample:")
        print("  python script.py trainval_paths_no-rd.txt standard.xyz ref.xyz")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_standard = sys.argv[2]
    output_ref = sys.argv[3]
    
    # Check input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    process_file_list(input_file, output_standard, output_ref)


if __name__ == "__main__":
    main()
