import os
import numpy as np
from ase.io import read, write
from pathlib import Path
import shutil

def voigt_6_to_full_3x3_stress(stress_vector):
    """
    Form a 3x3 stress matrix from a 6 component vector in Voigt notation
    """
    s1, s2, s3, s4, s5, s6 = stress_vector
    return np.array([[s1, s6, s5],
                     [s6, s2, s4],
                     [s5, s4, s3]])


def process_atoms(atoms):
    """
    Process a single Atoms object to convert fields from REF_ format to standard format
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    
    # Get the info dictionary
    info = atoms.info.copy()
    
    # 1. Convert REF_forces to forces in arrays
    forces = None
    if 'REF_forces' in atoms.arrays:
        forces = atoms.arrays['REF_forces'].copy()
        del atoms.arrays['REF_forces']
    
    # Store values we need
    energy_val = None
    free_energy_val = None
    stress_val = None
    
    # 2. Get energy from atomization_energy
    if 'atomization_energy' in info:
        energy_val = info['atomization_energy']
        del info['atomization_energy']
    
    # 3. Remove REF_energy
    if 'REF_energy' in info:
        del info['REF_energy']
    
    # Get free_energy if it exists. I don't think this happens
    if 'free_energy' in info:
        free_energy_val = info['free_energy']
        del info['free_energy']
    
    # 4. Convert REF_stress from Voigt to 3x3 matrix
    stress_voigt = None
    if 'REF_stress' in info:
        stress_voigt = info['REF_stress']
        ## Handle both string and array formats
        #if isinstance(stress_voigt, str):
        #    stress_voigt = np.array([float(x) for x in stress_voigt.split()])
        #elif not isinstance(stress_voigt, np.ndarray):
        #    stress_voigt = np.array(stress_voigt)
        #
        ## Convert to 3x3
        #stress_3x3 = voigt_6_to_full_3x3_stress(stress_voigt)
        ## Flatten to 9 components for storage
        #stress_val = stress_3x3.flatten()
        del info['REF_stress']
    
    # Update the atoms info with REF_stress (which should be written to info dict)
    #if stress_val is not None:
    #    info['REF_stress'] = stress_val
    
    new_info = {'energy' : energy_val,
                'free_energy': atoms.calc.results['free_energy']} 
    atoms.info = new_info
    
    # Set up calculator with energy and free_energy
    # This ensures they're stored properly in calculator results
    calc_results = {}
    #if energy_val is not None:
    #    calc_results['energy'] = energy_val
    #if free_energy_val is not None:
    #    calc_results['free_energy'] = atoms.calc.results['free_energy']
    if stress_voigt is not None:
        calc_results['stress'] = stress_voigt
    calc_results['forces'] = forces
    
    if calc_results:
        calc = SinglePointCalculator(atoms, **calc_results)
        atoms.calc = calc
    
    return atoms


#def process_atoms(atoms):
#    """
#    Process a single Atoms object to convert fields from REF_ format to standard format
#    """
#    # Get the info dictionary
#    info = atoms.info.copy()
#    
#    # 1. Convert REF_forces to forces in arrays
#    if 'REF_forces' in atoms.arrays:
#        atoms.arrays['forces'] = atoms.arrays['REF_forces'].copy()
#        del atoms.arrays['REF_forces']
#    
#    # Store values we need to reorder
#    energy_val = None
#    free_energy_val = None
#    stress_val = None
#    
#    # 2. Replace atomization_energy with energy
#    if 'atomization_energy' in info:
#        energy_val = info['atomization_energy']
#        del info['atomization_energy']
#    
#    # 3. Remove REF_energy
#    if 'REF_energy' in info:
#        del info['REF_energy']
#    
#    # Store free_energy if it exists
#    if 'free_energy' in info:
#        free_energy_val = info['free_energy']
#        del info['free_energy']
#    
#    # 4. Convert REF_stress from Voigt to 3x3 matrix
#    if 'REF_stress' in info:
#        stress_voigt = info['REF_stress']
#        # Handle both string and array formats
#        if isinstance(stress_voigt, str):
#            stress_voigt = np.array([float(x) for x in stress_voigt.split()])
#        elif not isinstance(stress_voigt, np.ndarray):
#            stress_voigt = np.array(stress_voigt)
#        
#        # Convert to 3x3
#        stress_3x3 = voigt_6_to_full_3x3_stress(stress_voigt)
#        # Flatten to 9 components for storage
#        stress_val = stress_3x3.flatten()
#        del info['REF_stress']
#    
#    # Now add them back in the desired order: energy, free_energy, REF_stress
#    new_info = {}
#    if energy_val is not None:
#        new_info['energy'] = energy_val
#    if free_energy_val is not None:
#        new_info['free_energy'] = free_energy_val
#    if stress_val is not None:
#        new_info['REF_stress'] = stress_val
#    
#    # Update the atoms info
#    atoms.info = new_info
#    
#    return atoms

#def process_atoms(atoms):
#    """
#    Process a single Atoms object to convert fields from REF_ format to standard format
#    """
#    # Get the info dictionary
#    info = atoms.info.copy()
#    
#    # 1. Convert REF_forces to forces in arrays
#    if 'REF_forces' in atoms.arrays:
#        atoms.arrays['forces'] = atoms.arrays['REF_forces'].copy()
#        del atoms.arrays['REF_forces']
#    
#    # 2. Replace atomization_energy with energy
#    if 'atomization_energy' in info:
#        info['energy'] = info['atomization_energy']
#        del info['atomization_energy']
#    
#    # 3. Remove REF_energy
#    if 'REF_energy' in info:
#        del info['REF_energy']
#    
#    # 4. Convert REF_stress from Voigt to 3x3 matrix
#    if 'REF_stress' in info:
#        #stress_voigt = info['REF_stress']
#        ## Handle both string and array formats
#        #if isinstance(stress_voigt, str):
#        #    stress_voigt = np.array([float(x) for x in stress_voigt.split()])
#        #elif not isinstance(stress_voigt, np.ndarray):
#        #    stress_voigt = np.array(stress_voigt)
#        #
#        ## Convert to 3x3
#        #stress_3x3 = voigt_6_to_full_3x3_stress(stress_voigt)
#        ## Flatten to 9 components for storage
#        #info['REF_stress'] = stress_3x3.flatten()
#        stress_voigt = info['REF_stress']
#        info['stress'] = stress_voigt 
#        del info['REF_stress']
#    
#    new_info = {'energy' : info['energy'],
#                'free_energy' : atoms.calc.results['free_energy'],
#                'stress' : info['stress']
#    }
#    # Update the atoms info
#    atoms.info = new_info
#    
#    return atoms

def process_xyz_file(input_path, output_path):
    """
    Process a single XYZ file and write the converted version
    """
    try:
        # Read all configurations from the file
        configs = read(input_path, index=':')
        
        # Handle single configuration case
        if not isinstance(configs, list):
            configs = [configs]
        
        # Process each configuration
        processed_configs = []
        for atoms in configs:
            processed_atoms = process_atoms(atoms)
            processed_configs.append(processed_atoms)
        
        # Write all configurations to output file
        write(output_path, processed_configs)
        print(f"Processed: {input_path} -> {output_path}")
        
    except Exception as e:
        print(f"Error processing {input_path}: {str(e)}")

def mirror_and_process_directory(source_dir, target_dir):
    """
    Mirror directory structure and process all XYZ files
    """
    source_path = Path(source_dir)
    target_path = Path(target_dir)
    
    # Create target directory if it doesn't exist
    target_path.mkdir(parents=True, exist_ok=True)
    
    # Walk through source directory
    for root, dirs, files in os.walk(source_path):
        # Get relative path from source
        rel_path = Path(root).relative_to(source_path)
        target_root = target_path / rel_path
        
        # Create subdirectories in target
        for dir_name in dirs:
            (target_root / dir_name).mkdir(parents=True, exist_ok=True)
        
        # Process XYZ files
        for file_name in files:
            if file_name.endswith('.xyz'):
                source_file = Path(root) / file_name
                target_file = target_root / file_name
                process_xyz_file(source_file, target_file)

def main():
    # Define source and target directories
    #source_dir = "processed_and_split_high-force-rm"
    #target_dir = "for_pod_split_high-force_rm"
    source_dir = "nothing_above_5"
    target_dir = "for_pod_nothing_above_5"

    
    print(f"Starting processing from '{source_dir}' to '{target_dir}'...")
    print("-" * 60)
    
    # Check if source directory exists
    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' does not exist!")
        return
    
    # Process all files
    mirror_and_process_directory(source_dir, target_dir)
    
    print("-" * 60)
    print("Processing complete!")

if __name__ == "__main__":
    main()
