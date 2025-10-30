#!/usr/bin/env python3
"""
Script to verify that configurations in processed_and_split files correspond
one-to-one with configurations in raw files by comparing free_energy and forces.

sonnet 4.0
"""

import os
from pathlib import Path
from ase.io import read
import argparse
import numpy as np
from collections import defaultdict

def extract_config_signature(atoms, is_raw=False):
    """
    Extract identifying signature from an ASE Atoms object.
    
    Args:
        atoms: ASE Atoms object
        is_raw: True if from raw files, False if from processed files
    
    Returns:
        tuple: (free_energy, forces_hash) for comparison
    """
    try:
        # Extract free energy
        free_energy = atoms.calc.results.get('free_energy', None)
        if free_energy is None:
            return None
        
        # Extract forces - different locations for raw vs processed
        if is_raw:
            # For raw files: forces are in atoms.calc.results or atoms.info
            if hasattr(atoms, 'calc') and atoms.calc is not None:
                forces = atoms.calc.results.get('forces', None)
            else:
                forces = atoms.arrays.get('forces', None)
        else:
            # For processed files: forces are in REF_forces
            forces = atoms.arrays.get('REF_forces', None)
        
        if forces is None:
            return None
            
        # Create a hash of forces for comparison (rounded to handle numerical precision)
        forces_rounded = np.round(forces, decimals=10)  # 10 decimal places precision
        forces_hash = hash(forces_rounded.tobytes())
        
        return (float(free_energy), forces_hash)
        
    except Exception as e:
        print(f"    Warning: Error extracting signature - {e}")
        return None

def load_config_signatures(file_path, is_raw=False):
    """
    Load configuration signatures from an XYZ file.
    
    Returns:
        list: List of (index, signature) tuples
    """
    print(f"  Loading configurations from {file_path}")
    
    try:
        configs = read(file_path, index=':')
        signatures = []
        
        for i, atoms in enumerate(configs):
            sig = extract_config_signature(atoms, is_raw=is_raw)
            if sig is not None:
                signatures.append((i, sig))
            else:
                print(f"    Warning: Could not extract signature for config {i}")
                
        print(f"    Loaded {len(signatures)} valid signatures out of {len(configs)} configs")
        return signatures
        
    except Exception as e:
        print(f"    Error loading {file_path}: {e}")
        return []

def check_file_correspondence(raw_file, trainval_file, test_file, tolerance=1e-8):
    """
    Check that configurations in trainval and test files correspond to raw file.
    """
    print(f"\nChecking correspondence:")
    print(f"  Raw: {raw_file}")
    print(f"  Trainval: {trainval_file}")
    print(f"  Test: {test_file}")
    
    # Load signatures from all files
    raw_sigs = load_config_signatures(raw_file, is_raw=True)
    trainval_sigs = load_config_signatures(trainval_file, is_raw=False) if trainval_file.exists() else []
    test_sigs = load_config_signatures(test_file, is_raw=False) if test_file.exists() else []
    
    if not raw_sigs:
        print("  ❌ No valid signatures found in raw file")
        return False
    
    # Create sets of signatures for comparison
    raw_sig_set = set(sig[1] for sig in raw_sigs)
    trainval_sig_set = set(sig[1] for sig in trainval_sigs)
    test_sig_set = set(sig[1] for sig in test_sigs)
    combined_processed_set = trainval_sig_set | test_sig_set
    
    print(f"  Raw configs: {len(raw_sig_set)}")
    print(f"  Trainval configs: {len(trainval_sig_set)}")
    print(f"  Test configs: {len(test_sig_set)}")
    print(f"  Total processed: {len(combined_processed_set)}")
    
    # Check for perfect correspondence
    missing_from_processed = raw_sig_set - combined_processed_set
    extra_in_processed = combined_processed_set - raw_sig_set
    overlap_trainval_test = trainval_sig_set & test_sig_set
    
    success = True
    
    if missing_from_processed:
        print(f"  ❌ {len(missing_from_processed)} configurations from raw file NOT found in processed files")
        success = False
        
    if extra_in_processed:
        print(f"  ❌ {len(extra_in_processed)} configurations in processed files NOT found in raw file")
        success = False
        
    if overlap_trainval_test:
        print(f"  ❌ {len(overlap_trainval_test)} configurations appear in BOTH trainval and test (data leakage!)")
        success = False
        
    if len(raw_sig_set) != len(combined_processed_set):
        print(f"  ❌ Count mismatch: {len(raw_sig_set)} raw vs {len(combined_processed_set)} processed")
        success = False
    
    if success:
        print("  ✅ Perfect correspondence - all raw configs found exactly once in processed files")
    
    return success

def find_corresponding_files(base_raw_dir, base_processed_dir, materials):
    """
    Find corresponding raw and processed file triplets.
    """
    file_triplets = []
    
    for material in materials:
        raw_material_dir = base_raw_dir / material
        processed_material_dir = base_processed_dir / material
        
        if not raw_material_dir.exists():
            print(f"Warning: Raw directory {raw_material_dir} does not exist")
            continue
            
        if not processed_material_dir.exists():
            print(f"Warning: Processed directory {processed_material_dir} does not exist")
            continue
        
        # Find all raw XYZ files
        raw_files = list(raw_material_dir.glob('*.xyz'))
        
        for raw_file in raw_files:
            # Construct expected processed filenames
            base_name = raw_file.stem
            trainval_file = processed_material_dir / 'trainval' / f"{base_name}_trainval.xyz"
            test_file = processed_material_dir / 'test' / f"{base_name}_test.xyz"
            
            # Check if both processed files exist
            if trainval_file.exists() and test_file.exists():
                file_triplets.append((raw_file, trainval_file, test_file))
            else:
                print(f"Warning: Missing processed files for {raw_file}")
                if not trainval_file.exists():
                    print(f"  Missing: {trainval_file}")
                if not test_file.exists():
                    print(f"  Missing: {test_file}")
    
    return file_triplets

def main():
    parser = argparse.ArgumentParser(description='Check correspondence between raw and processed configurations')
    parser.add_argument('--raw-dir', type=str, default='./raw',
                       help='Base directory containing raw files (default: ./raw)')
    parser.add_argument('--processed-dir', type=str, default='./processed_and_split',
                       help='Base directory containing processed files (default: ./processed_and_split)')
    parser.add_argument('--materials', nargs='+', default=['Hf', 'HfOx'],
                       help='Material folders to check (default: Hf HfOx)')
    parser.add_argument('--tolerance', type=float, default=1e-8,
                       help='Numerical tolerance for comparisons (default: 1e-8)')
    
    args = parser.parse_args()
    
    base_raw_dir = Path(args.raw_dir)
    base_processed_dir = Path(args.processed_dir)
    
    print(f"Checking configuration correspondence")
    print(f"Raw directory: {base_raw_dir}")
    print(f"Processed directory: {base_processed_dir}")
    print(f"Materials: {args.materials}")
    
    if not base_raw_dir.exists():
        print(f"Error: Raw directory {base_raw_dir} does not exist")
        return
        
    if not base_processed_dir.exists():
        print(f"Error: Processed directory {base_processed_dir} does not exist")
        return
    
    # Find all corresponding file triplets
    file_triplets = find_corresponding_files(base_raw_dir, base_processed_dir, args.materials)
    
    if not file_triplets:
        print("No corresponding file triplets found!")
        return
    
    print(f"\nFound {len(file_triplets)} file triplets to check")
    
    # Check each triplet
    total_checked = 0
    total_successful = 0
    
    for raw_file, trainval_file, test_file in file_triplets:
        success = check_file_correspondence(raw_file, trainval_file, test_file, args.tolerance)
        total_checked += 1
        if success:
            total_successful += 1
    
    print(f"\n{'='*80}")
    print(f"FINAL SUMMARY")
    print(f"{'='*80}")
    print(f"Total file triplets checked: {total_checked}")
    print(f"Successful correspondences: {total_successful}")
    print(f"Failed correspondences: {total_checked - total_successful}")
    
    if total_successful == total_checked:
        print("🎉 All files show perfect correspondence!")
    else:
        print("⚠️  Some files have correspondence issues - check output above")
    
    print("\nNote: This script compares free_energy and forces to verify correspondence.")
    print("Make sure your raw files have 'forces' and processed files have 'REF_forces'.")

if __name__ == "__main__":
    main()
