#!/usr/bin/env python3
"""
Sobol Sensitivity Analysis for holes_plate using Monte Carlo 
Alternative Version using pyapprox's sampling_based_sobol_indices function

This script uses pyapprox's AynchModel for PARALLEL execution of an external
simulation code  and computes Sobol sensitivity
indices using the built-in sampling_based_sobol_indices function from pyapprox.

Configuration:
- 50 total samples (executed in parallel)
- Monte Carlo method: pyapprox sampling_based_sobol_indices (Saltelli scheme)
- MAX_EVAL_CONCURRENCY: 2 (PARALLEL with 2 threads ONLY)

Monte Carlo Details:
- Uses pyapprox's sampling_based_sobol_indices function
- Sobol quasirandom sampling for better convergence
- Automatically handles Saltelli scheme internally
- Computes main effects and total effects

Usage:
    python sobol_parallel_mc.py

Requirements:
    - pyapprox installed
    - holes_plate executable
    - run_holes_plate.sh script
    - holeplate.template and hole.summit files
"""

import os
import sys
import numpy as np
from scipy import stats

from pyapprox.variables import IndependentMarginalsVariable
from pyapprox.interface.async_model import AynchModel

# Import the sampling_based_sobol_indices function from pyapprox
# First add the reference implementation to the path
ref_pyapprox_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "ref_pyapprox")
)
sys.path.insert(0, ref_pyapprox_path)

from pyapprox.analysis.sensitivity_analysis import (
    sampling_based_sobol_indices,
    get_isotropic_anova_indices,
)

# =============================================================================
# Configuration
# =============================================================================

# Path to shell script that runs holes_plate
SHELL_COMMAND = "./run_holes_plate.sh"

# Files that need to be linked into each work directory
LINK_FILENAMES = [
    os.path.abspath("holeplate.template"),
    os.path.abspath("hole.summit"),
    os.path.abspath("run_holes_plate.sh"),
]

# Parallel execution settings
MAX_EVAL_CONCURRENCY = 10  

# Monte Carlo settings
NUM_MC_SAMPLES = 1000  # Base samples for Monte Carlo

# Output directory for work directories
WORKDIR_BASENAME = "./workdirs/eval"

# =============================================================================
# Define Uncertain Variables
# =============================================================================


def create_variable():
    """
    Define the uncertain input variables and their distributions.

    Returns
    -------
    variable : IndependentMarginalsVariable
        Joint distribution of independent input parameters
    var_names : list
        Names of variables for labeling results
    """
    # Uniform distributions for uncertain parameters
    univariate_variables = [
        stats.uniform(loc=4.36e3, scale=(4.84e3-4.36e3) ),  # rho: [0.8, 1.2]
        stats.uniform(loc=90e9, scale=(138e9-90e9) ),  # E: [1.5e11, 2.5e11] Pa
        stats.uniform(loc=0.35, scale=0.37-0.35),  # nu: [0.25, 0.35]
    ]

    variable = IndependentMarginalsVariable(univariate_variables)
    var_names = ["rho", "E", "nu"]

    return variable, var_names


# =============================================================================
# Create Model Wrapper
# =============================================================================


def create_asynch_model(workdir_basename=None):
    """
    Create an AynchModel wrapper for parallel holes_plate evaluation.

    Parameters
    ----------
    workdir_basename : str, optional
        Base path for work directories. If None, uses temp directories.

    Returns
    -------
    model : AynchModel
        Callable model that evaluates holes_plate in parallel
    """
    if workdir_basename is not None:
        # Ensure parent directory exists
        parent_dir = os.path.dirname(workdir_basename)
        if parent_dir and not os.path.exists(parent_dir):
            os.makedirs(parent_dir)

    model = AynchModel(
        shell_command=SHELL_COMMAND,
        max_eval_concurrency=MAX_EVAL_CONCURRENCY,
        workdir_basename=workdir_basename,
        link_filenames=LINK_FILENAMES,
        params_filename="params.in",
        results_filename="results.out",
        save_workdirs="limited",  # Save only params/results, delete rest
        model_name="holes_plate",
    )

    return model


# =============================================================================
# Monte Carlo Sobol Analysis using pyapprox function
# =============================================================================


def run_sobol_analysis():
    """
    Run Sobol sensitivity analysis for holes_plate using pyapprox's
    sampling_based_sobol_indices function.

    Returns
    -------
    results : dict
        Dictionary containing sensitivity analysis results
    """
    print("=" * 70)
    print("Sobol Sensitivity Analysis - Parallel Monte Carlo (Version 2)")
    print("Using pyapprox's sampling_based_sobol_indices function")
    print("=" * 70)

    # 1. Create variable definitions
    print("\n[1/4] Defining uncertain variables...")
    variable, var_names = create_variable()
    print(f"  Variables: {var_names}")
    print(f"  Distributions: Uniform")
    nvars = len(var_names)

    # 2. Create model wrapper
    print("\n[2/4] Creating parallel model wrapper...")
    print(
        f"  Max concurrent evaluations: {MAX_EVAL_CONCURRENCY}"
    )
    model = create_asynch_model(workdir_basename=WORKDIR_BASENAME)

    # 3. Define interaction terms
    # get_isotropic_anova_indices(nvars, order) generates all interaction terms
    # up to the specified order. Order 2 includes main effects and pairwise interactions.
    print(f"\n[3/4] Defining interaction terms...")
    max_order = 2  # Main effects and pairwise interactions
    interaction_terms = get_isotropic_anova_indices(nvars, max_order)
    nterms = interaction_terms.shape[1]
    print(f"  Max interaction order: {max_order}")
    print(f"  Number of interaction terms: {nterms}")
    print(f"  Interaction terms shape: {interaction_terms.shape}")

    # Print interaction terms for reference
    print(f"\n  Interaction terms (each column is one term):")
    for i in range(nterms):
        active_vars = np.where(interaction_terms[:, i] > 0)[0]
        if len(active_vars) == 1:
            print(f"    Term {i}: Variable {active_vars[0]} (main effect)")
        else:
            print(f"    Term {i}: Variables {list(active_vars)} (interaction)")

    # 4. Run sampling-based Sobol analysis
    print(f"\n[4/4] Running sampling-based Sobol analysis...")
    print(f"  Base samples: {NUM_MC_SAMPLES}")
    print(f"  Sampling method: Sobol quasirandom sequence")

    # Calculate total number of model evaluations
    # sampling_based_sobol_indices uses:
    # - 2 sample sets A and B (2 * nsamples)
    # - nterms additional sample sets (nterms * nsamples)
    # Total = (2 + nterms) * nsamples
    total_evals = (2 + nterms) * NUM_MC_SAMPLES
    print(f"  Total model evaluations: {total_evals}")
    print(f"    = (2 base sets + {nterms} interaction sets) × {NUM_MC_SAMPLES} samples")
    print(f"\n  Starting parallel evaluations...")

    # Call pyapprox's sampling_based_sobol_indices function
    sobol_indices, total_effects, variance, mean = sampling_based_sobol_indices(
        fun=model,
        variables=variable,
        interaction_terms=interaction_terms,
        nsamples=NUM_MC_SAMPLES,
        sampling_method="sobol",  # Use Sobol quasirandom sequence
        qmc_start_index=0,
    )

    print("  Sobol indices computed successfully")

    # Extract main effects (first-order indices)
    # Main effects are the terms where only one variable is active
    main_effect_indices = []
    for i in range(nterms):
        if interaction_terms[:, i].sum() == 1:
            main_effect_indices.append(i)

    main_effects = sobol_indices[main_effect_indices, :]

    # sobol_indices has shape (nterms, nqoi)
    # total_effects has shape (nvars, nqoi)
    num_qoi = sobol_indices.shape[1]

    # ==========================================================================
    # Display Results
    # ==========================================================================
    print("\n" + "=" * 70)
    print("RESULTS - Parallel Monte Carlo V2 (2 threads, 50 samples)")
    print("Using pyapprox's sampling_based_sobol_indices function")
    print("=" * 70)

    # Assuming QoI names from the template
    qoi_names = ["E_effective", "stress", "strain", "max_von_Mises_stress"]
    if len(qoi_names) != num_qoi:
        qoi_names = [f"QoI_{i}" for i in range(num_qoi)]

    print("\n--- Main Effects (First-order Sobol Indices) ---")
    print(f"{'Variable':<10}", end="")
    for qoi in qoi_names:
        print(f"{qoi:<20}", end="")
    print()

    for i, var_name in enumerate(var_names):
        print(f"{var_name:<10}", end="")
        for j in range(num_qoi):
            print(f"{main_effects[i, j]:<20.8f}", end="")
        print()

    print("\n--- Total Effects ---")
    print(f"{'Variable':<10}", end="")
    for qoi in qoi_names:
        print(f"{qoi:<20}", end="")
    print()

    for i, var_name in enumerate(var_names):
        print(f"{var_name:<10}", end="")
        for j in range(num_qoi):
            print(f"{total_effects[i, j]:<20.8f}", end="")
        print()

    # Display all Sobol indices (including interactions)
    print("\n--- All Sobol Indices (including interactions) ---")
    for i in range(nterms):
        active_vars = np.where(interaction_terms[:, i] > 0)[0]
        var_label = ",".join([var_names[v] for v in active_vars])
        print(f"\nTerm {i}: ({var_label})")
        print(f"{'QoI':<20}", end="")
        print(f"{'Sobol Index':<15}")
        for j in range(num_qoi):
            print(f"{qoi_names[j]:<20}", end="")
            print(f"{sobol_indices[i, j]:<20.8f}")

    # Display mean and variance
    print("\n--- Mean and Variance ---")
    print(f"{'QoI':<20}", end="")
    print(f"{'Mean':<20}", end="")
    print(f"{'Variance':<20}")
    for j in range(num_qoi):
        print(f"{qoi_names[j]:<20}", end="")
        print(f"{mean[j]:<20.8e}", end="")
        print(f"{variance[j]:<20.8e}")

    # Return results for further analysis
    results = {
        "variable": variable,
        "var_names": var_names,
        "qoi_names": qoi_names,
        "interaction_terms": interaction_terms,
        "sobol_indices": sobol_indices,
        "main_effects": main_effects,
        "total_effects": total_effects,
        "variance": variance,
        "mean": mean,
    }

    # Save results
    output_file = "sobol_parallel_mc_results.npz"
    np.savez(
        output_file,
        interaction_terms=interaction_terms,
        sobol_indices=sobol_indices,
        main_effects=main_effects,
        total_effects=total_effects,
        variance=variance,
        mean=mean,
    )
    print(f"\nResults saved to {output_file}")

    return results


# =============================================================================
# Entry Point # =============================================================================
if __name__ == "__main__":
    # Ensure OMP_NUM_THREADS=1 for best parallel performance
    os.environ["OMP_NUM_THREADS"] = "1"

    print("\nStarting Sobol Sensitivity Analysis (Parallel Monte Carlo V2)")
    print(f"  Configuration: 50 samples, pyapprox function, parallel with 2 threads\n")

    results = run_sobol_analysis()

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)
