#!/usr/bin/env python3
"""
Split an XYZ file into training and validation sets based on index files.
"""

from ase.io import read, write

def read_indices(filename):
    """Read indices from a text file (one per line)."""
    with open(filename, 'r') as f:
        indices = [int(line.strip()) for line in f if line.strip()]
    return indices

def main():
    # Read all configurations from the XYZ file
    print("Reading configurations from trainval.xyz...")
    configs = read('nothing_above_5/all/trainval/trainval.xyz', index=':')
    print(f"Total configurations: {len(configs)}")
    
    # Read training indices
    print("\nReading training indices...")
    train_indices = read_indices('new_dpp_train_indices.txt')
    print(f"Training indices: {len(train_indices)}")
    
    # Read validation indices
    print("Reading validation indices...")
    val_indices = read_indices('new_dpp_validation_indices.txt')
    print(f"Validation indices: {len(val_indices)}")
    
    # Extract training configurations
    print("\nExtracting training configurations...")
    train_configs = [configs[i] for i in train_indices]
    
    # Extract validation configurations
    print("Extracting validation configurations...")
    val_configs = [configs[i] for i in val_indices]
    
    # Write training set
    print(f"\nWriting {len(train_configs)} configurations to new_dpp_train.xyz...")
    write('new_dpp_train.xyz', train_configs)
    
    # Write validation set
    print(f"Writing {len(val_configs)} configurations to new_dpp_validation.xyz...")
    write('new_dpp_validation.xyz', val_configs)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
