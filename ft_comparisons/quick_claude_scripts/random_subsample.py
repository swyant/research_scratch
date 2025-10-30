import os
import random
from pathlib import Path
from ase.io import read, write

# Set random seed for reproducibility (optional)
random.seed(42)

# Define input directories
dirs = [
    'processed_and_split_high-force-rm/Hf/trainval',
    'processed_and_split_high-force-rm/HfOx/trainval'
]

# Output file
output_file = './random_subsample_trainval-1.xyz'

# Sampling fraction
sample_fraction = 0.1

# Clear output file if it exists
#if os.path.exists(output_file):
#    os.remove(output_file)

# Process each directory
for directory in dirs:
    dir_path = Path(directory)
    
    # Find all xyz files in directory
    xyz_files = list(dir_path.glob('*.xyz'))
    
    print(f"\nProcessing {directory}:")
    print(f"Found {len(xyz_files)} xyz files")
    
    # Process each xyz file
    for xyz_file in xyz_files:
        # Read all configurations from file
        atoms_list = read(xyz_file, index=':')
        
        # Calculate number of configurations to sample
        n_total = len(atoms_list)
        n_sample = max(1, int(n_total * sample_fraction))
        
        # Randomly sample configurations
        sampled_atoms = random.sample(atoms_list, n_sample)
        
        # Append to output file
        write(output_file, sampled_atoms, append=True)
        
        print(f"  {xyz_file.name}: sampled {n_sample}/{n_total} configurations")

print(f"\nAll sampled configurations written to {output_file}")
