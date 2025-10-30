#!/usr/bin/env python3
"""
Sanity check script for train/test splits.
Verifies that atomization energies (or other properties) don't appear in both
trainval and test sets, which would indicate data leakage.
"""

import os
from pathlib import Path
from ase.io import read
import argparse
from collections import defaultdict
import numpy as np

def extract_properties(atoms_list, property_keys=None):
    """
    Extract specified properties from a list of ASE Atoms objects.
    
    Args:
        atoms_list: List of ASE Atoms objects
        property_keys: List of property keys to extract. If None, tries common ones.
    
    Returns:
        dict: Dictionary mapping property names to lists of values
    """
    if property_keys is None:
        # Common property keys to check
        property_keys = [
            'atomization_energy',
            'energy', 
            'total_energy',
            'formation_energy',
            'DFT_energy',
            'energy_per_atom'
        ]
    
    properties = defaultdict(list)
    
    for atoms in atoms_list:
        info = atoms.info
        for key in property_keys:
            if key in info:
                properties[key].append(info[key])
    
    return dict(properties)

def check_overlap(trainval_props, test_props, tolerance=1e-10):
    """
    Check for overlapping values between trainval and test property sets.
    
    Args:
        trainval_props: Dict of property arrays from trainval set
        test_props: Dict of property arrays from test set  
        tolerance: Numerical tolerance for floating point comparison
    
    Returns:
        dict: Results of overlap analysis
    """
    results = {}
    
    for prop_name in trainval_props:
        if prop_name not in test_props:
            results[prop_name] = {
                'overlap_count': 0,
                'message': f'Property {prop_name} not found in test set'
            }
            continue
            
        trainval_values = np.array(trainval_props[prop_name])
        test_values = np.array(test_props[prop_name])
        
        # Find overlapping values within tolerance
        overlaps = []
        for i, test_val in enumerate(test_values):
            matches = np.abs(trainval_values - test_val) <= tolerance
            if np.any(matches):
                overlaps.append((i, test_val, trainval_values[matches]))
        
        results[prop_name] = {
            'trainval_count': len(trainval_values),
            'test_count': len(test_values),
            'overlap_count': len(overlaps),
            'overlaps': overlaps[:10] if overlaps else [],  # Show first 10
            'message': f'Found {len(overlaps)} overlapping values' if overlaps else 'No overlaps found'
        }
    
    return results

def analyze_file_pair(trainval_file, test_file, tolerance=1e-10):
    """
    Analyze a pair of trainval/test files for overlapping properties.
    """
    print(f"\nAnalyzing file pair:")
    print(f"  Trainval: {trainval_file}")
    print(f"  Test: {test_file}")
    
    try:
        # Read configurations
        trainval_configs = read(trainval_file, index=':')
        test_configs = read(test_file, index=':')
        
        print(f"  Loaded {len(trainval_configs)} trainval configs, {len(test_configs)} test configs")
        
        # Extract properties
        trainval_props = extract_properties(trainval_configs)
        test_props = extract_properties(test_configs)
        
        if not trainval_props and not test_props:
            print("  Warning: No recognized properties found in either file")
            return
        
        # Check for overlaps
        results = check_overlap(trainval_props, test_props, tolerance)
        
        # Report results
        any_overlaps = False
        for prop_name, result in results.items():
            print(f"  {prop_name}: {result['message']}")
            if result['overlap_count'] > 0:
                any_overlaps = True
                print(f"    ⚠️  WARNING: {result['overlap_count']} overlapping values found!")
                if result['overlaps']:
                    print(f"    First few overlaps: {[ov[1] for ov in result['overlaps'][:3]]}")
        
        if not any_overlaps:
            print("  ✅ No overlaps detected - split looks good!")
            
    except Exception as e:
        print(f"  ❌ Error analyzing files: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Sanity check train/test splits for data leakage')
    #parser.add_argument('--base-dir', type=str, default='./processed_and_split',
    #                   help='Base directory containing split files (default: ./processed_and_split)')
    parser.add_argument('--base-dir', type=str, default='processed_and_split_high-force-rm',
                       help='Base directory containing split files (default: processed_and_split_high-force-rm)')
    parser.add_argument('--tolerance', type=float, default=1e-10,
                       help='Numerical tolerance for property comparison (default: 1e-10)')
    parser.add_argument('--materials', nargs='+', default=['Hf', 'HfOx'],
                       help='Material folders to check (default: Hf HfOx)')
    
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    
    if not base_dir.exists():
        print(f"Error: Base directory {base_dir} does not exist")
        return
    
    print(f"Checking train/test splits in {base_dir}")
    print(f"Numerical tolerance: {args.tolerance}")
    
    total_files_checked = 0
    problems_found = 0
    
    for material in args.materials:
        material_dir = base_dir / material
        trainval_dir = material_dir / 'trainval'
        test_dir = material_dir / 'test'
        
        print(f"\n{'='*60}")
        print(f"CHECKING MATERIAL: {material}")
        print(f"{'='*60}")
        
        if not trainval_dir.exists():
            print(f"Warning: {trainval_dir} does not exist")
            continue
        if not test_dir.exists():
            print(f"Warning: {test_dir} does not exist")
            continue
        
        # Find matching trainval/test file pairs
        trainval_files = list(trainval_dir.glob('*_trainval.xyz'))
        
        if not trainval_files:
            print(f"No trainval files found in {trainval_dir}")
            continue
        
        for trainval_file in sorted(trainval_files):
            # Construct corresponding test filename
            base_name = trainval_file.stem.replace('_trainval', '')
            test_filename = f"{base_name}_test.xyz"
            test_file = test_dir / test_filename
            
            if not test_file.exists():
                print(f"Warning: Expected test file {test_file} not found")
                continue
            
            # Analyze the file pair
            analyze_file_pair(trainval_file, test_file, args.tolerance)
            total_files_checked += 1
    
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Total file pairs checked: {total_files_checked}")
    
    if total_files_checked == 0:
        print("⚠️  No file pairs found to check!")
    else:
        print("✅ Sanity check complete!")
        print("\nIf any overlaps were found above, you may have data leakage.")
        print("Consider checking your splitting logic or using a different random seed.")

if __name__ == "__main__":
    main()
