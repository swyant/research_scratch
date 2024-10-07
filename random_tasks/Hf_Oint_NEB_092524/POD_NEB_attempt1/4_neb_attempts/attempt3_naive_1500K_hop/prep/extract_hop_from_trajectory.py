#!/usr/bin/env python 
import ase.io 
from ase import Atoms


spec_map = {1: "Hf", 2: "O"}

traj_file = "../../../../../comb_HF-O_diffusivity_072924/redo_process_with_POD/4_run_NVT/1500K/1500K_nvt.lammpstrj"

full_traj = ase.io.read(traj_file, format="lammps-dump-text", index=":")

for sys in full_traj:
    sys.set_chemical_symbols([spec_map[number] for number in sys.numbers])

jump_images = full_traj[1498:1509]

for i, image in enumerate(jump_images):
    ase.io.write(f"1500K_image.coords.{i}", image, format="lammps-data", specorder=["Hf", "O"], masses=True)
