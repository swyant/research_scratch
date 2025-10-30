#!/usr/bin/env python3
"""
Script to process XYZ files from ./raw/Hf and ./raw/HfOx folders.
Transforms energy and forces fields according to specifications.

Initial draft with sonnet 4.0
"""

import os
import glob
from pathlib import Path
from ase.io import read, write
#from ase import units
from ase.units import create_units
from ase.calculators.singlepoint import SinglePointCalculator

"""
Dionysios: The energies are: Hf = -112.28902966 Ry (sigma = 0.00035 Ryd) and O = -32.54698947 Ry (sigma = 0.000005 Ryd)
"""

# Reference atomic energies in Rydberg

# The values Dio gave me
#E0_Hf_Ry = -112.28902966
#E0_O_Ry = -32.54698947

# The values I back-computed with 2022 CODATA values
#E0_Hf_Ry = -112.41606520152305989309
#E0_O_Ry = -32.62034915207989080779

# According to Dio, I need to use the 2006 CODATA values. 
# Truly doubt this affects anything
units_2006 = create_units('2006')

# The values I back-computed with 2006 CODATA values 
E0_Hf_Ry = -112.41607411999959316330 
E0_O_Ry = -32.6203517400001970321 

# Convert to eV using ASE units
E0_Hf_eV = E0_Hf_Ry * units_2006['Ry']  # Rydberg to eV conversion
E0_O_eV = E0_O_Ry * units_2006['Ry']

print(f"Reference energies:")
print(f"E0_Hf: {E0_Hf_eV:.8f} eV")
print(f"E0_O: {E0_O_eV:.8f} eV")

def process_atoms_configuration(atoms, config_index):
    """
    Process a single atoms configuration.
    
    Args:
        atoms: ASE Atoms object
        config_index: Index of configuration for logging
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Check if we have a calculator with results
        if not hasattr(atoms, 'calc') or atoms.calc is None:
            print(f"    Config {config_index}: No calculator found")
            return False
            
        if not hasattr(atoms.calc, 'results') or not atoms.calc.results:
            print(f"    Config {config_index}: No calculator results found")
            return False
            
        calc_results = atoms.calc.results
        
        # Get original energy from calculator
        if 'energy' not in calc_results:
            print(f"    Config {config_index}: No 'energy' in calculator results")
            return False
            
        original_energy = calc_results['energy']
        
        # Count atoms for REF_energy calculation
        n_Hf = sum(1 for symbol in atoms.get_chemical_symbols() if symbol == 'Hf')
        n_O = sum(1 for symbol in atoms.get_chemical_symbols() if symbol == 'O')
        
        # Calculate REF_energy by adding back isolated atomic energies
        ref_energy = original_energy + n_Hf * E0_Hf_eV + n_O * E0_O_eV
        
        #print(f"    Config {config_index}: Atoms: {n_Hf} Hf, {n_O} O")
        #print(f"    Config {config_index}: Original energy: {original_energy:.8f} eV")
        #print(f"    Config {config_index}: REF_energy: {ref_energy:.8f} eV")
        
        # 1. Store original energy as atomization_energy in atoms.info
        atoms.info['atomization_energy'] = original_energy
        
        # 2. Store REF_energy in atoms.info
        atoms.info['REF_energy'] = ref_energy
        
        # 3. Move forces from calculator to atoms.arrays as REF_forces
        if 'forces' in calc_results:
            atoms.arrays['REF_forces'] = calc_results['forces'].copy()
            #print(f"    Config {config_index}: Moved forces to REF_forces")
        
        # 4. Move stress from calculator to atoms.info as REF_stress
        if 'stress' in calc_results:
            atoms.info['REF_stress'] = calc_results['stress'].copy()
            #print(f"    Config {config_index}: Moved stress to REF_stress")
        
        # 5. Create new calculator with only free_energy (and any other fields we want to keep)
        new_results = {}
        
        # Keep free_energy if it exists
        if 'free_energy' in calc_results:
            print(f"free energy - reference_energy = {calc_results['free_energy']-ref_energy:.8f}")
            new_results['free_energy'] = calc_results['free_energy']
            #print(f"    Config {config_index}: Keeping free_energy: {calc_results['free_energy']:.8f} eV")
        
        # Keep any other fields that aren't energy, forces, or stress
        excluded_keys = {'energy', 'forces', 'stress'}
        for key, value in calc_results.items():
            if key not in excluded_keys and key != 'free_energy':  # free_energy already handled
                new_results[key] = value
                #print(f"    Config {config_index}: Keeping calculator field: {key}")
        
        # Create new SinglePoint calculator with filtered results
        if new_results:
            new_calc = SinglePointCalculator(atoms, **new_results)
            atoms.calc = new_calc
        else:
            # If no results to keep, remove calculator entirely
            atoms.calc = None
            #print(f"    Config {config_index}: Removed calculator (no fields to keep)")
        
        return True
        
    except Exception as e:
        print(f"    Config {config_index}: Error processing: {str(e)}")
        return False

def process_xyz_file(input_path, output_path):
    """
    Process a single XYZ file (which may contain multiple configurations).
    
    Args:
        input_path (str): Path to input XYZ file
        output_path (str): Path to output XYZ file
    """
    try:
        # Read all configurations from the file
        all_atoms = read(input_path, index=':')  # ':' reads all configurations
        
        # Handle case where there's only one configuration (returns single Atoms object)
        if not isinstance(all_atoms, list):
            all_atoms = [all_atoms]
        
        print(f"Processing {input_path}:")
        print(f"  Found {len(all_atoms)} configurations")
        
        processed_atoms = []
        successful_configs = 0
        
        # Process each configuration
        for i, atoms in enumerate(all_atoms):
            if process_atoms_configuration(atoms, i+1):
                processed_atoms.append(atoms)
                successful_configs += 1
            else:
                print(f"    Config {i+1}: Skipping due to processing error")
        
        if not processed_atoms:
            print(f"  No configurations could be processed successfully")
            return False
        
        #print(f"  Successfully processed {successful_configs}/{len(all_atoms)} configurations")
        
        # Write all processed configurations to output file
        write(output_path, processed_atoms, format='extxyz')
        #print(f"  Saved to: {output_path}")
        
        return True
        
    except Exception as e:
        print(f"Error processing {input_path}: {str(e)}")
        return False


def main():
    """Main processing function."""
    
    # Define input and output directories
    input_dirs = ['./raw/Hf', './raw/HfOx']
    output_dirs = ['./processed/Hf', './processed/HfOx']
    
    # Create output directories if they don't exist
    for output_dir in output_dirs:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    total_processed = 0
    total_files = 0
    
    # Process each directory
    for input_dir, output_dir in zip(input_dirs, output_dirs):
        if not os.path.exists(input_dir):
            print(f"Warning: Input directory {input_dir} does not exist")
            continue
            
        print(f"\nProcessing directory: {input_dir}")
        print(f"Output directory: {output_dir}")
        
        # Find all top-level xyz files
        xyz_pattern = os.path.join(input_dir, '*.xyz')
        xyz_files = glob.glob(xyz_pattern)
        
        if not xyz_files:
            print(f"No XYZ files found in {input_dir}")
            continue
            
        print(f"Found {len(xyz_files)} XYZ files")
        
        # Process each file
        for xyz_file in xyz_files:
            total_files += 1
            filename = os.path.basename(xyz_file)
            output_path = os.path.join(output_dir, filename)
            
            if process_xyz_file(xyz_file, output_path):
                total_processed += 1
            
            print()  # Add blank line for readability
    
    # Summary
    print("="*50)
    print(f"Processing complete!")
    print(f"Total files found: {total_files}")
    print(f"Successfully processed: {total_processed}")
    print(f"Failed: {total_files - total_processed}")

if __name__ == "__main__":
    main()
