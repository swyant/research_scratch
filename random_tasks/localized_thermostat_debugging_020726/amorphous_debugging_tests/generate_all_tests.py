"""Generate all amorphous debugging test cases.

Run with: uv run python amorphous_debugging_tests/generate_all_tests.py
"""

import os
import sys
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from localized_thermostat.config import ThermostatConfig
from localized_thermostat.run_simulation import generate_files, EXAMPLE_DIR, DEFAULT_DATA_FILE, POTENTIAL_FILES

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def copy_common_files(output_dir):
    """Copy DATA and POD files to output directory."""
    os.makedirs(output_dir, exist_ok=True)
    # Copy DATA file
    data_src = DEFAULT_DATA_FILE
    data_dst = os.path.join(output_dir, "thermostat.data")
    if not os.path.exists(data_dst):
        shutil.copy2(data_src, data_dst)
    # Copy POD files
    for fname in POTENTIAL_FILES:
        src = os.path.join(EXAMPLE_DIR, fname)
        dst = os.path.join(output_dir, fname)
        if not os.path.exists(dst):
            shutil.copy2(src, dst)


def write_script(output_dir, script_content):
    """Write LAMMPS input script."""
    path = os.path.join(output_dir, "thermostat.in")
    with open(path, 'w') as f:
        f.write(script_content)
    print(f"  Written: {path}")


def full_dump_commands(every, diag_dir="diagnostics"):
    """Return dump commands that output BOTH mobile-only and full-structure trajectories."""
    return f"""\
# --- Diagnostics ---
shell mkdir {diag_dir}
compute ke_per_atom all ke/atom

# Mobile atoms only (backward compatible)
dump diag_traj mobile custom {every} {diag_dir}/ke_trajectory.lammpstrj id type x y z c_ke_per_atom
dump_modify diag_traj sort id format line "%d %d %14.8f %14.8f %14.8f %14.8e"

# FULL structure including frozen atoms (for visualization)
dump full_traj all custom {every} {diag_dir}/full_trajectory.lammpstrj id type x y z c_ke_per_atom
dump_modify full_traj sort id format line "%d %d %14.8f %14.8f %14.8f %14.8e"
"""


# ============================================================
# Test 01: Bare stability — NVE (no thermostat, no frozen atoms)
# ============================================================
def generate_01_bare_nve():
    name = "01_bare_stability_nve"
    output_dir = os.path.join(BASE_DIR, name)
    copy_common_files(output_dir)
    os.makedirs(os.path.join(output_dir, "diagnostics"), exist_ok=True)

    script = f"""\
# Bare stability test: NVE, no thermostat, no frozen atoms
# Tests whether the POD potential + structure are stable on their own
# All atoms initialized at 50 K

units metal
atom_style atomic
boundary p p p

read_data thermostat.data

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

pair_style pod
pair_coeff * * HfO2_param.pod HfO2_coefficients.pod Hf O

# Initialize all atoms at low temperature
velocity all create 50.0 12345 dist gaussian

# Pure NVE — no thermostat
fix integ all nve
timestep 0.001

{full_dump_commands(100)}

# Thermo
compute temp_all all temp
thermo 100
thermo_style custom step temp c_temp_all etotal pe ke press vol

run 5000
"""
    write_script(output_dir, script)
    print(f"  [{name}] Bare NVE at 50K, 5000 steps, all atoms mobile")


# ============================================================
# Test 02: Bare stability — NVT at 50K (global Langevin, no frozen)
# ============================================================
def generate_02_bare_nvt():
    name = "02_bare_stability_nvt_50K"
    output_dir = os.path.join(BASE_DIR, name)
    copy_common_files(output_dir)
    os.makedirs(os.path.join(output_dir, "diagnostics"), exist_ok=True)

    script = f"""\
# Bare stability test: NVT with global Langevin at 50K
# Tests whether structure stays crystalline under standard thermostatting
# No frozen atoms, no localized thermostat

units metal
atom_style atomic
boundary p p p

read_data thermostat.data

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

pair_style pod
pair_coeff * * HfO2_param.pod HfO2_coefficients.pod Hf O

# Initialize all atoms at 50 K
velocity all create 50.0 12345 dist gaussian

# NVE + Langevin = NVT
fix integ all nve
fix lang all langevin 50.0 50.0 0.1 54321
timestep 0.001

{full_dump_commands(100)}

# Thermo
compute temp_all all temp
thermo 100
thermo_style custom step temp c_temp_all etotal pe ke press vol

run 5000
"""
    write_script(output_dir, script)
    print(f"  [{name}] Bare NVT at 50K, 5000 steps, all atoms mobile")


# ============================================================
# Test 03: Single shell, mode 1, 50K
# ============================================================
def generate_03_single_shell_mode1():
    name = "03_single_shell_mode1_50K"
    output_dir = os.path.join(BASE_DIR, name)

    config = ThermostatConfig(
        core_radius=5.0,
        transition_width=10.0,
        n_shells=1,           # Single shell
        core_temp=50.0,
        edge_temp=50.0,       # Same temp everywhere = uniform 50K
        mode=1,
        timestep=0.001,
        thermostat_damp=0.1,
        n_steps=5000,
        diag_every=100,
        output_dir=output_dir,
    )

    generate_files(config)
    # Now patch the script to add full-structure dump
    _patch_add_full_dump(output_dir)
    print(f"  [{name}] Mode 1, 1 shell, uniform 50K, frozen atoms present")


# ============================================================
# Test 04: Single shell, mode 3, 50K
# ============================================================
def generate_04_single_shell_mode3():
    name = "04_single_shell_mode3_50K"
    output_dir = os.path.join(BASE_DIR, name)

    config = ThermostatConfig(
        core_radius=5.0,
        transition_width=10.0,
        n_shells=1,
        core_temp=50.0,
        edge_temp=50.0,
        mode=3,
        timestep=0.001,
        thermostat_damp=0.1,
        n_steps=5000,
        diag_every=100,
        output_dir=output_dir,
    )

    generate_files(config)
    _patch_add_full_dump(output_dir)
    print(f"  [{name}] Mode 3, 1 shell, uniform 50K, frozen atoms present")


# ============================================================
# Test 05: Uniform thermostat 50K (all mobile, no frozen, with Langevin per-zone)
# This isolates the frozen boundary effect — if this is crystalline but
# test 03 is amorphous, the frozen atoms are the problem.
# ============================================================
def generate_05_uniform_thermostat():
    name = "05_uniform_thermostat_50K"
    output_dir = os.path.join(BASE_DIR, name)
    copy_common_files(output_dir)
    os.makedirs(os.path.join(output_dir, "diagnostics"), exist_ok=True)

    # Use localized thermostat but with a huge radius so nothing is frozen
    # core_radius + transition_width = 50 Å > half box → everything is mobile
    script = f"""\
# Uniform thermostat at 50K: localized thermostat with no frozen atoms
# core_radius=50 → entire box is "core" → all atoms thermostated at 50K
# This tests whether the localized thermostat machinery causes issues
# independent of the frozen boundary

units metal
atom_style atomic
boundary p p p

read_data thermostat.data

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

pair_style pod
pair_coeff * * HfO2_param.pod HfO2_coefficients.pod Hf O

# All atoms are mobile (no frozen group)
group mobile id *

# Initialize all at 50K
velocity all create 50.0 12345 dist gaussian

fix integ mobile nve
timestep 0.001

fix lang_all mobile langevin 50.0 50.0 0.1 54321

{full_dump_commands(100)}

compute temp_all all temp
thermo 100
thermo_style custom step temp c_temp_all etotal pe ke press vol

run 5000
"""
    write_script(output_dir, script)
    print(f"  [{name}] Uniform 50K thermostat, no frozen atoms")


# ============================================================
# Test 06: Full 4-shell thermostat, mode 1, LOW temperatures (50K core, 10K edge)
# ============================================================
def generate_06_full_thermostat_low_temp():
    name = "06_full_thermostat_mode1_low_temp"
    output_dir = os.path.join(BASE_DIR, name)

    config = ThermostatConfig(
        core_radius=5.0,
        transition_width=10.0,
        n_shells=4,
        core_temp=50.0,       # Much lower than the 500K used before
        edge_temp=10.0,
        mode=1,
        timestep=0.001,
        thermostat_damp=0.1,
        n_steps=5000,
        diag_every=100,
        output_dir=output_dir,
    )

    generate_files(config)
    _patch_add_full_dump(output_dir)
    print(f"  [{name}] Mode 1, 4 shells, 50K core / 10K edge")


# ============================================================
# Test 07: Timestep sensitivity — same as test 03 but with 0.0005 ps timestep
# If the structure is unstable with 0.001 ps but stable with 0.0005 ps,
# the timestep is too large for POD forces.
# ============================================================
def generate_07_timestep_sensitivity():
    name = "07_timestep_sensitivity"
    output_dir = os.path.join(BASE_DIR, name)

    config = ThermostatConfig(
        core_radius=5.0,
        transition_width=10.0,
        n_shells=1,
        core_temp=50.0,
        edge_temp=50.0,
        mode=1,
        timestep=0.0005,      # Half the default timestep
        thermostat_damp=0.1,
        n_steps=10000,        # Double steps to cover same simulation time
        diag_every=200,
        output_dir=output_dir,
    )

    generate_files(config)
    _patch_add_full_dump(output_dir)
    print(f"  [{name}] Mode 1, 1 shell, 50K, timestep=0.0005 ps (half normal)")


def _patch_add_full_dump(output_dir):
    """Add full-structure dump command to an existing thermostat.in file."""
    script_path = os.path.join(output_dir, "thermostat.in")
    with open(script_path) as f:
        content = f.read()

    # Add full trajectory dump after the mobile-only dump
    full_dump_line = (
        "\n# FULL structure including frozen atoms (for visualization)\n"
        "dump full_traj all custom 100 diagnostics/full_trajectory.lammpstrj "
        "id type x y z c_ke_per_atom\n"
        'dump_modify full_traj sort id format line "%d %d %14.8f %14.8f %14.8f %14.8e"\n'
    )

    # Insert after the existing dump_modify line
    marker = 'dump_modify diag_traj'
    if marker in content:
        idx = content.index(marker)
        end_of_line = content.index('\n', idx)
        content = content[:end_of_line+1] + full_dump_line + content[end_of_line+1:]

    with open(script_path, 'w') as f:
        f.write(content)


if __name__ == '__main__':
    print("Generating amorphous debugging test cases...\n")

    generate_01_bare_nve()
    generate_02_bare_nvt()
    generate_03_single_shell_mode1()
    generate_04_single_shell_mode3()
    generate_05_uniform_thermostat()
    generate_06_full_thermostat_low_temp()
    generate_07_timestep_sensitivity()

    print(f"\nAll test cases generated in {BASE_DIR}/")
    print("\nTo run a test:")
    print("  cd amorphous_debugging_tests/<test_dir>")
    print("  /opt/lammps/build/lmp -in thermostat.in -log log.lammps")
    print("\nTo analyze results:")
    print("  uv run python amorphous_debugging_tests/analyze_crystallinity.py \\")
    print("    <test_dir>/diagnostics/full_trajectory.lammpstrj \\")
    print("    --data-file example_simulation/pod_relaxed_Hf_supercell_Oint_nonvt_DATA \\")
    print("    --all-frames")
