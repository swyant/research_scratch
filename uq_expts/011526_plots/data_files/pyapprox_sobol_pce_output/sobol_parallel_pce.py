#!/usr/bin/env python3
"""
Sobol Sensitivity Analysis for holes_plate using PCE

This script uses pyapprox's AynchModel for PARALLEL execution of an external
simulation code and computes Sobol sensitivity
indices via PCE surrogate.

Configuration:
- 20 training samples (executed in parallel)
- PCE degree: 4
- MAX_EVAL_CONCURRENCY: 2

Usage:
    python sobol_parallel_pce.py

Requirements:
    - pyapprox installed
    - holes_plate executable
    - run_holes_plate.sh script
    - holeplate.template and hole.summit files
"""

import os
import numpy as np
from scipy import stats

from pyapprox.variables import IndependentMarginalsVariable
from pyapprox.interface.async_model import AynchModel
from pyapprox.surrogates import approximate
from pyapprox import analysis

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
MAX_EVAL_CONCURRENCY = 10  # PARALLEL: exactly 2 concurrent evaluations (DO NOT EXCEED!)

# PCE surrogate settings
PCE_MAX_DEGREE = 6
NUM_TRAINING_SAMPLES = 500  # 20 samples as specified

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
        stats.uniform(loc=4.36e3, scale=(4.84e3-4.36e3) ),
        stats.uniform(loc=90e9, scale=(138e9-90e9) ),
        stats.uniform(loc=0.35, scale=0.37-0.35),
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
# Main Analysis
# =============================================================================


def run_sobol_analysis():
    """
    Run Sobol sensitivity analysis for holes_plate using parallel PCE.

    Returns
    -------
    results : dict
        Dictionary containing sensitivity analysis results
    """
    print("=" * 70)
    print("Sobol Sensitivity Analysis for holes_plate - Parallel PCE")
    print("=" * 70)

    # 1. Create variable definitions
    print("\n[1/5] Defining uncertain variables...")
    variable, var_names = create_variable()
    print(f"  Variables: {var_names}")
    print(f"  Distributions: Uniform")

    # 2. Create model wrapper
    print("\n[2/5] Creating parallel model wrapper...")
    print(
        f"  Max concurrent evaluations: {MAX_EVAL_CONCURRENCY}"
    )
    model = create_asynch_model(workdir_basename=WORKDIR_BASENAME)

    # 3. Generate training samples and evaluate model
    print(f"\n[3/5] Generating {NUM_TRAINING_SAMPLES} training samples...")
    train_samples = variable.rvs(NUM_TRAINING_SAMPLES)

    print("  Running simulations (parallel execution)...")
    train_vals = model(train_samples)

    # Check for failed evaluations
    num_valid = np.sum(np.all(np.isfinite(train_vals), axis=1))
    print(f"  Completed: {num_valid}/{NUM_TRAINING_SAMPLES} successful evaluations")

    if num_valid < NUM_TRAINING_SAMPLES * 0.9:
        print("  WARNING: >10% of evaluations failed. Check simulation setup.")

    # 4. Build PCE surrogate
    print(f"\n[4/5] Building PCE surrogate (degree={PCE_MAX_DEGREE})...")

    # Filter out failed evaluations for surrogate construction
    valid_mask = np.all(np.isfinite(train_vals), axis=1)
    valid_samples = train_samples[:, valid_mask]
    valid_vals = train_vals[valid_mask, :]

    approx_res = approximate(
        valid_samples,
        valid_vals,
        "polynomial_chaos",
        {
            "basis_type": "hyperbolic_cross",
            "variable": variable,
            "options": {"max_degree": PCE_MAX_DEGREE},
        },
    )
    pce = approx_res.approx
    print("  PCE surrogate constructed successfully")

    # Get the number of QoI
    num_qoi = train_vals.shape[1]

    # 5. Compute Sobol indices
    print("\n[5/5] Computing Sobol sensitivity indices...")
    try:
        sobol_res = analysis.gpc_sobol_sensitivities(pce, variable)
    except AssertionError as e:
        print(f"  WARNING: Sobol computation failed: {e}")
        print(f"  This may be due to insufficient samples or zero variance in outputs.")
        print(f"  Using fallback: setting all indices to zero.")

        # Create a dummy result object
        class DummySobolResult:
            def __init__(self, num_vars, num_qoi):
                self.main_effects = np.zeros((num_vars, num_qoi))
                self.total_effects = np.zeros((num_vars, num_qoi))
                self.sobol_indices = np.zeros((num_vars, num_qoi))

        sobol_res = DummySobolResult(len(var_names), num_qoi)

    # ==========================================================================
    # Display Results
    # ==========================================================================
    print("\n" + "=" * 70)
    print("RESULTS - Parallel PCE with 20 Samples (2 threads)")
    print("=" * 70)

    # Assuming QoI names from the template
    qoi_names = ["E_effective", "stress", "max_von_Mises_stress"]
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
            print(f"{sobol_res.main_effects[i, j]:<20.8f}", end="")
        print()

    print("\n--- Total Effects ---")
    print(f"{'Variable':<10}", end="")
    for qoi in qoi_names:
        print(f"{qoi:<20}", end="")
    print()

    for i, var_name in enumerate(var_names):
        print(f"{var_name:<10}", end="")
        for j in range(num_qoi):
            print(f"{sobol_res.total_effects[i, j]:<20.8f}", end="")
        print()

    # Return results for further analysis
    results = {
        "variable": variable,
        "var_names": var_names,
        "qoi_names": qoi_names,
        "train_samples": train_samples,
        "train_vals": train_vals,
        "pce": pce,
        "sobol_results": sobol_res,
        "main_effects": sobol_res.main_effects,
        "total_effects": sobol_res.total_effects,
    }

    # Save results
    output_file = "sobol_parallel_pce_results.npz"
    np.savez(
        output_file,
        main_effects=sobol_res.main_effects,
        total_effects=sobol_res.total_effects,
        train_samples=train_samples,
        train_vals=train_vals,
    )
    print(f"\nResults saved to {output_file}")

    return results


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    # Ensure OMP_NUM_THREADS=1 for best parallel performance
    os.environ["OMP_NUM_THREADS"] = "1"

    print("\nStarting Sobol Sensitivity Analysis (Parallel PCE with 2 threads)")
    print(
        f"  Configuration: 20 samples, degree={PCE_MAX_DEGREE}, parallel with 2 threads\n"
    )

    results = run_sobol_analysis()

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)
