#!/usr/bin/env python3
"""
Script to process extended XYZ files with multiple configurations:
1. Randomize configuration order
2. Split into 85% train/val and 15% test sets
3. Save to organized folder structure

Sonnet 4.0 9/25/25

pre-filtering force w/ Sonnet 4.5
"""

import os
import random
from pathlib import Path
from ase.io import read, write
import argparse
import numpy as np

def check_force_threshold(config, threshold=15.0):
    """
    Check if any force component exceeds the threshold.
    
    Args:
        config: ASE Atoms object
        threshold: Maximum allowed absolute force value in eV/Å
    
    Returns:
        bool: True if forces are within threshold, False if any exceed
    """
    if 'REF_forces' not in config.arrays:
        return True  # No forces present, include config
    
    forces = config.arrays['REF_forces']
    max_force = np.abs(forces).max()
    
    return max_force <= threshold


def process_xyz_file(input_file, output_base_dir, material_folder, force_threshold=15.0):
    """
    Process a single XYZ file: randomize, split, and save configurations.
    
    Args:
        input_file (Path): Path to the input XYZ file
        output_base_dir (Path): Base output directory (./processed_and_split)
        material_folder (str): Material folder name (Hf or HfOx)
        force_threshold (float): Maximum allowed absolute force value in eV/Å
    """
    print(f"Processing {input_file}...")
    
    try:
        # Read all configurations from the XYZ file
        configurations = read(input_file, index=':')
        total_configs = len(configurations)
        print(f"  Found {total_configs} configurations")        

        if len(configurations) == 0:
            print(f"  Warning: No configurations found in {input_file}")
            return

        # Filter configurations based on force threshold
        filtered_configs = []
        high_force_configs = []
        
        for config in configurations:
            if check_force_threshold(config, force_threshold):
                filtered_configs.append(config)
            else:
                high_force_configs.append(config)
        
        num_filtered = len(filtered_configs)
        num_high_force = len(high_force_configs)
        
        print(f"  Force filtering results:")
        print(f"    - Kept: {num_filtered} configs (within ±{force_threshold} eV/Å)")
        print(f"    - Removed: {num_high_force} configs (exceeding threshold)")

        if num_filtered == 0:
            print(f"  Warning: No configurations passed force filtering!")
            # Still save high-force configs if they exist
            if num_high_force > 0:
                material_output_dir = output_base_dir / material_folder
                high_force_dir = material_output_dir / 'high_force'
                high_force_dir.mkdir(parents=True, exist_ok=True)

                stem = input_file.stem
                suffix = input_file.suffix
                high_force_filename = f"{stem}_high-force{suffix}"
                high_force_path = high_force_dir / high_force_filename

                write(high_force_path, high_force_configs)
                print(f"  Wrote {num_high_force} high-force configs to {high_force_path}")
            return

        # Randomize the order of configurations
        random.shuffle(filtered_configs)
        
        # Calculate split sizes
        train_val_size = int(0.85 * num_filtered)
        test_size = num_filtered - train_val_size
        
        print(f"  Split: {train_val_size} train/val, {test_size} test")
        
        # Split configurations
        train_val_configs = filtered_configs[:train_val_size]
        test_configs = filtered_configs[train_val_size:]
        
        # Prepare output paths
        material_output_dir = output_base_dir / material_folder
        trainval_dir = material_output_dir / 'trainval'
        test_dir = material_output_dir / 'test'
        high_force_dir = material_output_dir / 'high_force'
        
        # Create directories if they don't exist
        trainval_dir.mkdir(parents=True, exist_ok=True)
        test_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate output filenames
        stem = input_file.stem  # filename without extension
        suffix = input_file.suffix  # .xyz
        
        trainval_filename = f"{stem}_trainval{suffix}"
        test_filename = f"{stem}_test{suffix}"
        
        trainval_path = trainval_dir / trainval_filename
        test_path = test_dir / test_filename
        
        # Write the split files
        if train_val_configs:
            write(trainval_path, train_val_configs)
            print(f"  Wrote {len(train_val_configs)} configs to {trainval_path}")
        
        if test_configs:
            write(test_path, test_configs)
            print(f"  Wrote {len(test_configs)} configs to {test_path}")

        # Write high-force configurations if any exist
        if high_force_configs:
            high_force_dir.mkdir(parents=True, exist_ok=True)
            high_force_filename = f"{stem}_high-force{suffix}"
            high_force_path = high_force_dir / high_force_filename

            write(high_force_path, high_force_configs)
            print(f"  Wrote {num_high_force} high-force configs to {high_force_path}")
            
    except Exception as e:
        print(f"  Error processing {input_file}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Process and split XYZ configuration files')
    parser.add_argument('--seed', type=int, default=42, 
                       help='Random seed for reproducible splits (default: 42)')
    parser.add_argument('--input-base', type=str, default='./processed',
                       help='Base input directory (default: ./processed)')
    parser.add_argument('--output-base', type=str, default='./processed_and_split_high-force-rm',
                       help='Base output directory (default: ./processed_and_split_high-force-rm)')
    parser.add_argument('--force-threshold', type=float, default=15.0,
                       help='Maximum allowed absolute force value in eV/Å (default: 15.0)')
    
    args = parser.parse_args()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    
    # Define paths
    input_base_dir = Path(args.input_base)
    output_base_dir = Path(args.output_base)
    
    # Check if input directories exist
    if not input_base_dir.exists():
        print(f"Error: Input directory {input_base_dir} does not exist")
        return

    print(f"Configuration:")
    print(f"  Force threshold: ±{args.force_threshold} eV/Å")
    print(f"  Random seed: {args.seed}")
    print(f"  Input directory: {input_base_dir}")
    print(f"  Output directory: {output_base_dir}")

    # Process both material folders
    material_folders = ['Hf', 'HfOx']

    for material in material_folders:
        input_material_dir = input_base_dir / material

        if not input_material_dir.exists():
            print(f"\nWarning: {input_material_dir} does not exist, skipping...")
            continue

        print(f"\n{'='*60}")
        print(f"Processing material folder: {material}")
        print(f"{'='*60}")
        print(f"Input directory: {input_material_dir}")

        # Find all .xyz files in the material directory
        xyz_files = list(input_material_dir.glob('*.xyz'))

        if not xyz_files:
            print(f"  No .xyz files found in {input_material_dir}")
            continue

        print(f"Found {len(xyz_files)} XYZ files")

        # Process each XYZ file
        for xyz_file in sorted(xyz_files):
            process_xyz_file(xyz_file, output_base_dir, material, args.force_threshold)

    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"{'='*60}")
    print(f"Output structure:")
    print(f"  {output_base_dir}/Hf/trainval/")
    print(f"  {output_base_dir}/Hf/test/")
    print(f"  {output_base_dir}/Hf/high_force/")
    print(f"  {output_base_dir}/HfOx/trainval/")
    print(f"  {output_base_dir}/HfOx/test/")
    print(f"  {output_base_dir}/HfOx/high_force/")


if __name__ == "__main__":
    main()
