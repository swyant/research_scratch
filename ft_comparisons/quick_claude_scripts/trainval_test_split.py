#!/usr/bin/env python3
"""
Script to process extended XYZ files with multiple configurations:
1. Randomize configuration order
2. Split into 85% train/val and 15% test sets
3. Save to organized folder structure

Sonnet 4.0 9/25/25
"""

import os
import random
from pathlib import Path
from ase.io import read, write
import argparse

def process_xyz_file(input_file, output_base_dir, material_folder):
    """
    Process a single XYZ file: randomize, split, and save configurations.
    
    Args:
        input_file (Path): Path to the input XYZ file
        output_base_dir (Path): Base output directory (./processed_and_split)
        material_folder (str): Material folder name (Hf or HfOx)
    """
    print(f"Processing {input_file}...")
    
    try:
        # Read all configurations from the XYZ file
        configurations = read(input_file, index=':')
        print(f"  Found {len(configurations)} configurations")
        
        if len(configurations) == 0:
            print(f"  Warning: No configurations found in {input_file}")
            return
        
        # Randomize the order of configurations
        random.shuffle(configurations)
        
        # Calculate split sizes
        total_configs = len(configurations)
        train_val_size = int(0.85 * total_configs)
        test_size = total_configs - train_val_size
        
        print(f"  Split: {train_val_size} train/val, {test_size} test")
        
        # Split configurations
        train_val_configs = configurations[:train_val_size]
        test_configs = configurations[train_val_size:]
        
        # Prepare output paths
        material_output_dir = output_base_dir / material_folder
        trainval_dir = material_output_dir / 'trainval'
        test_dir = material_output_dir / 'test'
        
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
            
    except Exception as e:
        print(f"  Error processing {input_file}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Process and split XYZ configuration files')
    parser.add_argument('--seed', type=int, default=42, 
                       help='Random seed for reproducible splits (default: 42)')
    parser.add_argument('--input-base', type=str, default='./processed',
                       help='Base input directory (default: ./processed)')
    parser.add_argument('--output-base', type=str, default='./processed_and_split',
                       help='Base output directory (default: ./processed_and_split)')
    
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
    
    # Process both material folders
    material_folders = ['Hf', 'HfOx']
    
    for material in material_folders:
        input_material_dir = input_base_dir / material
        
        if not input_material_dir.exists():
            print(f"Warning: {input_material_dir} does not exist, skipping...")
            continue
            
        print(f"\nProcessing material folder: {material}")
        print(f"Input directory: {input_material_dir}")
        
        # Find all .xyz files in the material directory
        xyz_files = list(input_material_dir.glob('*.xyz'))
        
        if not xyz_files:
            print(f"  No .xyz files found in {input_material_dir}")
            continue
            
        print(f"  Found {len(xyz_files)} XYZ files")
        
        # Process each XYZ file
        for xyz_file in sorted(xyz_files):
            process_xyz_file(xyz_file, output_base_dir, material)
    
    print(f"\nProcessing complete!")
    print(f"Output structure:")
    print(f"  {output_base_dir}/Hf/trainval/")
    print(f"  {output_base_dir}/Hf/test/")
    print(f"  {output_base_dir}/HfOx/trainval/") 
    print(f"  {output_base_dir}/HfOx/test/")

if __name__ == "__main__":
    main()
