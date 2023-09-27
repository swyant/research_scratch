import ase.io

input_fname = "train.xyz"
configs = ase.io.read(input_fname, index=":")
for idx, config in enumerate(configs):
    ase.io.write(f"data_files/prl_siviraman_HfO2_DATA_{idx}", config, format="lammps-data", masses=True)
