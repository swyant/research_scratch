import ase.io

input_fname = "trial_dft_configs.xyz"
configs = ase.io.read(input_fname, index=":")
for idx, config in enumerate(configs):
    ase.io.write(f"data_files/lammps_pod_stress_test_DATA_{idx}", config, format="lammps-data", masses=True)
