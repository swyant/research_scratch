"""
Utility functions for evaluating MACE potentials using ASE.

This module provides functions to:
1. Load and label atomic configurations from directory structure
2. Compute energies and forces using calculators (e.g., MACE)
3. Calculate mean absolute errors for different groups/subgroups
4. Display results in a clean format

Sonnet 4.0
Note the error calculation is pretty redundant but for now this is fine
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Union
import numpy as np
from ase import Atoms
from ase.io import read
from ase.calculators.calculator import Calculator
import pandas as pd
from collections import defaultdict
import warnings


def load_configurations(data_dir: Union[str, Path]) -> Dict[str, List[Atoms]]:
    """
    Load all xyz files from the processed_and_split directory structure.
    
    Args:
        data_dir: Path to the "processed_and_split" directory
        
    Returns:
        Dictionary with keys as labels and values as lists of ASE Atoms objects.
        Labels format: "{material}_{split}_{temp}" where:
        - material: "Hf" or "HfOx"
        - split: "trainval" or "test"  
        - temp: "0K" or "no0K"
    """
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileNotFoundError(f"Data directory {data_path} does not exist")
    
    configurations = defaultdict(list)
    
    # Expected materials and splits
    materials = ["Hf", "HfOx"]
    splits = ["trainval", "test"]
    
    for material in materials:
        material_dir = data_path / material
        if not material_dir.exists():
            warnings.warn(f"Material directory {material_dir} not found, skipping")
            continue
            
        for split in splits:
            split_dir = material_dir / split
            if not split_dir.exists():
                warnings.warn(f"Split directory {split_dir} not found, skipping")
                continue
            
            # Find all xyz files in this directory
            xyz_files = list(split_dir.glob("*.xyz"))
            if not xyz_files:
                warnings.warn(f"No xyz files found in {split_dir}")
                continue
            
            print(f"Found {len(xyz_files)} files in {material}/{split}")
            
            for xyz_file in xyz_files:
                try:
                    # Read all configurations from the file
                    atoms_list = read(str(xyz_file), index=":")
                    if not isinstance(atoms_list, list):
                        atoms_list = [atoms_list]
                    
                    # Determine temperature category from filename
                    temp_category = "0K" if "0K" in xyz_file.name else "no0K"
                    
                    # Create label
                    label = f"{material}_{split}_{temp_category}"
                    
                    # Add configurations to the dictionary
                    configurations[label].extend(atoms_list)
                    
                    print(f"  {xyz_file.name}: {len(atoms_list)} configs -> {label}")
                    
                except Exception as e:
                    warnings.warn(f"Error reading {xyz_file}: {e}")
                    continue
    
    # Convert defaultdict to regular dict and print summary
    configurations = dict(configurations)
    print("\nConfiguration summary:")
    total_configs = 0
    for label, atoms_list in configurations.items():
        print(f"  {label}: {len(atoms_list)} configurations")
        total_configs += len(atoms_list)
    print(f"Total configurations loaded: {total_configs}")
    
    return configurations


def compute_predictions(configurations: Dict[str, List[Atoms]], 
                       calculator: Calculator) -> Dict[str, List[Atoms]]:
    """
    Compute energies and forces for all configurations using the provided calculator.
    
    Args:
        configurations: Dictionary of labeled configuration lists from load_configurations
        calculator: ASE calculator (e.g., MACE calculator)
        
    Returns:
        Dictionary with same structure but containing copies of atoms with computed
        energies and forces from the calculator
    """
    predictions = {}
    
    total_configs = sum(len(atoms_list) for atoms_list in configurations.values())
    processed_configs = 0
    
    print(f"Computing predictions for {total_configs} configurations...")
    
    for label, atoms_list in configurations.items():
        print(f"\nProcessing {label}: {len(atoms_list)} configurations")
        pred_atoms_list = []
        
        for i, atoms in enumerate(atoms_list):
            try:
                # Create a copy to avoid modifying original
                atoms_copy = atoms.copy()
                
                # Set calculator and compute properties
                atoms_copy.calc = calculator
                
                # Force calculation of energy and forces
                energy = atoms_copy.get_potential_energy()
                forces = atoms_copy.get_forces()
                
                pred_atoms_list.append(atoms_copy)
                processed_configs += 1
                
                # Progress indicator
                if (i + 1) % 100 == 0 or (i + 1) == len(atoms_list):
                    print(f"  Processed {i + 1}/{len(atoms_list)} configs")
                    
            except Exception as e:
                warnings.warn(f"Error computing properties for config {i} in {label}: {e}")
                continue
        
        predictions[label] = pred_atoms_list
        print(f"  Successfully processed {len(pred_atoms_list)}/{len(atoms_list)} configs")
    
    print(f"\nTotal predictions computed: {processed_configs}/{total_configs}")
    return predictions


def _safe_get_energy(atoms: Atoms, is_reference: bool = False) -> Optional[float]:
    """Safely extract energy from atoms object."""
    try:
        if is_reference:
            # For reference data, look for REF_energy in info dict
            if 'REF_energy' in atoms.info:
                energy = atoms.info['REF_energy']
                # Normalize by number of atoms
                return energy / len(atoms)
            # Also check arrays in case it's stored there
            if hasattr(atoms, 'arrays') and 'REF_energy' in atoms.arrays:
                energy = atoms.arrays['REF_energy'][0]  # Assuming scalar
                return energy / len(atoms)
        else:
            # For predicted data, use standard methods
            if atoms.calc is not None and hasattr(atoms.calc, 'results'):
                if 'energy' in atoms.calc.results:
                    energy = atoms.calc.results['energy']
                    return energy / len(atoms)
            # Fallback methods
            if hasattr(atoms, 'get_potential_energy'):
                energy = atoms.get_potential_energy()
                return energy / len(atoms)
        return None
    except:
        return None


def _safe_get_forces(atoms: Atoms, is_reference: bool = False) -> Optional[np.ndarray]:
    """Safely extract forces from atoms object."""
    try:
        if is_reference:
            # For reference data, look for REF_forces in arrays
            if hasattr(atoms, 'arrays') and 'REF_forces' in atoms.arrays:
                return atoms.arrays['REF_forces']
            # Also check info dict
            if 'REF_forces' in atoms.info:
                return atoms.info['REF_forces']
        else:
            # For predicted data, use standard methods
            if atoms.calc is not None and hasattr(atoms.calc, 'results'):
                if 'forces' in atoms.calc.results:
                    return atoms.calc.results['forces']
            # Fallback methods
            if hasattr(atoms, 'get_forces'):
                return atoms.get_forces()
        return None
    except:
        return None


def calculate_errors(reference_configs: Dict[str, List[Atoms]], 
                    predicted_configs: Dict[str, List[Atoms]]) -> Dict[str, Dict[str, float]]:
    """
    Calculate mean absolute errors for energies and forces across different groups.
    
    Args:
        reference_configs: Reference configurations with DFT energies/forces
        predicted_configs: Predicted configurations from calculator
        
    Returns:
        Dictionary with error metrics for each group/subgroup:
        {
            "group_name": {
                "energy_mae": float,
                "force_mae": float,
                "n_configs": int,
                "n_atoms_total": int
            }
        }
    """
    # First collect errors for individual labels
    label_errors = {}
    all_energy_errors = []
    all_force_errors = []
    
    for label in reference_configs.keys():
        if label not in predicted_configs:
            warnings.warn(f"Label {label} not found in predictions, skipping")
            continue
            
        ref_atoms_list = reference_configs[label]
        pred_atoms_list = predicted_configs[label]
        
        if len(ref_atoms_list) != len(pred_atoms_list):
            warnings.warn(f"Mismatch in number of configs for {label}: "
                         f"{len(ref_atoms_list)} ref vs {len(pred_atoms_list)} pred")
            min_len = min(len(ref_atoms_list), len(pred_atoms_list))
            ref_atoms_list = ref_atoms_list[:min_len]
            pred_atoms_list = pred_atoms_list[:min_len]
        
        energy_errors = []
        force_errors = []
        n_atoms_total = 0
        
        for ref_atoms, pred_atoms in zip(ref_atoms_list, pred_atoms_list):
            # Energy errors (per atom)
            ref_energy = _safe_get_energy(ref_atoms, is_reference=True)
            pred_energy = _safe_get_energy(pred_atoms, is_reference=False)
            
            if ref_energy is not None and pred_energy is not None:
                energy_errors.append(abs(ref_energy - pred_energy))
                all_energy_errors.append(abs(ref_energy - pred_energy))
            
            # Force errors
            ref_forces = _safe_get_forces(ref_atoms, is_reference=True)
            pred_forces = _safe_get_forces(pred_atoms, is_reference=False)
            
            if ref_forces is not None and pred_forces is not None:
                if ref_forces.shape == pred_forces.shape:
                    # Flatten forces to get all components (fx, fy, fz for all atoms)
                    ref_forces_flat = ref_forces.flatten()
                    pred_forces_flat = pred_forces.flatten()
                    force_errors.extend(np.abs(ref_forces_flat - pred_forces_flat))
                    all_force_errors.extend(np.abs(ref_forces_flat - pred_forces_flat))
                    n_atoms_total += len(ref_atoms)
        
        # Store errors for this label
        label_errors[label] = {
            'energy_mae': np.mean(energy_errors) if energy_errors else np.nan,
            'force_mae': np.mean(force_errors) if force_errors else np.nan,
            'n_configs': len(ref_atoms_list),
            'n_atoms_total': n_atoms_total
        }
    
    # Now aggregate errors for different groups
    error_results = {}
    
    # Copy individual label results
    error_results.update(label_errors)
    
    # Define aggregation groups
    aggregation_groups = {
        'total_trainval': [k for k in label_errors.keys() if 'trainval' in k],
        'total_test': [k for k in label_errors.keys() if 'test' in k],
        'Hf_trainval': [k for k in label_errors.keys() if k.startswith('Hf_') and 'trainval' in k],
        'HfOx_trainval': [k for k in label_errors.keys() if k.startswith('HfOx_') and 'trainval' in k],
        'Hf_test': [k for k in label_errors.keys() if k.startswith('Hf_') and 'test' in k],
        'HfOx_test': [k for k in label_errors.keys() if k.startswith('HfOx_') and 'test' in k],
    }
    
    # Calculate aggregated errors
    for group_name, label_list in aggregation_groups.items():
        if not label_list:
            continue
            
        # Collect all errors for this group
        group_energy_errors = []
        group_force_errors = []
        group_n_configs = 0
        group_n_atoms = 0
        
        for label in label_list:
            if label not in reference_configs or label not in predicted_configs:
                continue
                
            ref_atoms_list = reference_configs[label]
            pred_atoms_list = predicted_configs[label]
            min_len = min(len(ref_atoms_list), len(pred_atoms_list))
            
            for i in range(min_len):
                ref_atoms = ref_atoms_list[i]
                pred_atoms = pred_atoms_list[i]
                
                # Energy
                ref_energy = _safe_get_energy(ref_atoms, is_reference=True)
                pred_energy = _safe_get_energy(pred_atoms, is_reference=False)
                if ref_energy is not None and pred_energy is not None:
                    group_energy_errors.append(abs(ref_energy - pred_energy))
                
                # Forces
                ref_forces = _safe_get_forces(ref_atoms, is_reference=True)
                pred_forces = _safe_get_forces(pred_atoms, is_reference=False)
                if ref_forces is not None and pred_forces is not None and ref_forces.shape == pred_forces.shape:
                    ref_forces_flat = ref_forces.flatten()
                    pred_forces_flat = pred_forces.flatten()
                    group_force_errors.extend(np.abs(ref_forces_flat - pred_forces_flat))
                    group_n_atoms += len(ref_atoms)
                    
                group_n_configs += 1
        
        error_results[group_name] = {
            'energy_mae': np.mean(group_energy_errors) if group_energy_errors else np.nan,
            'force_mae': np.mean(group_force_errors) if group_force_errors else np.nan,
            'n_configs': group_n_configs,
            'n_atoms_total': group_n_atoms
        }
    
    return error_results


def display_errors(error_results: Dict[str, Dict[str, float]], 
                  title: str = "Model Evaluation Results") -> None:
    """
    Display error results in a clean, formatted table.
    
    Args:
        error_results: Dictionary from calculate_errors function
        title: Title for the results table
    """
    print(f"\n{'=' * 60}")
    print(f"{title:^60}")
    print(f"{'=' * 60}")
    
    # Convert to DataFrame for nice formatting
    df_data = []
    for group_name, metrics in error_results.items():
        df_data.append({
            'Group': group_name,
            'Energy MAE (eV/atom)': f"{metrics['energy_mae']:.4f}" if not np.isnan(metrics['energy_mae']) else "N/A",
            'Force MAE (eV/Å)': f"{metrics['force_mae']:.4f}" if not np.isnan(metrics['force_mae']) else "N/A",
            'N Configs': metrics['n_configs'],
            'N Atoms': metrics['n_atoms_total']
        })
    
    df = pd.DataFrame(df_data)
    
    # Sort by logical order
    order = ['total_trainval', 'total_test', 
             'Hf_trainval', 'Hf_test', 
             'HfOx_trainval', 'HfOx_test',
             'Hf_trainval_no0K', 'Hf_trainval_0K',
             'Hf_test_no0K', 'Hf_test_0K',
             'HfOx_trainval_no0K', 'HfOx_trainval_0K',
             'HfOx_test_no0K', 'HfOx_test_0K']
    
    # Reorder DataFrame according to the specified order
    ordered_rows = []
    for group in order:
        mask = df['Group'] == group
        if mask.any():
            ordered_rows.append(df[mask])
    
    # Add any remaining rows not in the order
    remaining_groups = set(df['Group']) - set(order)
    for group in sorted(remaining_groups):
        mask = df['Group'] == group
        if mask.any():
            ordered_rows.append(df[mask])
    
    if ordered_rows:
        df_ordered = pd.concat(ordered_rows, ignore_index=True)
    else:
        df_ordered = df
    
    # Display with nice formatting
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    
    print(df_ordered.to_string(index=False))
    print(f"{'=' * 60}\n")


def evaluate_model(data_dir: Union[str, Path], 
                  calculator: Calculator,
                  title: str = "Model Evaluation Results") -> Dict[str, Dict[str, float]]:
    """
    Complete pipeline to evaluate a model.
    
    Args:
        data_dir: Path to processed_and_split directory
        calculator: ASE calculator to evaluate
        title: Title for results display
        
    Returns:
        Dictionary with error results
    """
    print("Step 1: Loading configurations...")
    reference_configs = load_configurations(data_dir)
    
    print("\nStep 2: Computing predictions...")
    predicted_configs = compute_predictions(reference_configs, calculator)
    
    print("\nStep 3: Calculating errors...")
    error_results = calculate_errors(reference_configs, predicted_configs)
    
    print("\nStep 4: Displaying results...")
    display_errors(error_results, title)
    
    return error_results


## Example usage function for testing
#def example_usage():
#    """
#    Example of how to use these utilities.
#    This function shows the typical workflow.
#    """
#    from mace.calculators import mace_mp
#    
#    # Example with MACE-MP calculator
#    calculator = mace_mp(model="medium", dispersion=False, default_dtype="float32")
#    
#    # Run evaluation
#    data_dir = "processed_and_split"  # Path to your data directory
#    results = evaluate_model(data_dir, calculator, "MACE-MP Medium Evaluation")
#    
#    return results


if __name__ == "__main__":
    # Run example if script is executed directly
    print("This is the utils.py module for MACE evaluation.")
    print("Import this module in your Jupyter notebook and use the functions:")
    print("  - load_configurations()")
    print("  - compute_predictions()")
    print("  - calculate_errors()")
    print("  - display_errors()")
    print("  - evaluate_model()  # Complete pipeline")
