#!/usr/bin/env python
from ase import Atoms
import ase.io 

relaxed_sys = ase.io.read("./hexag_Hf_0xz0yz_DATA", format="lammps-data")
relaxed_cell = relaxed_sys.cell

for i in range(7):
    infile  = f"OT_path_DATA{i}"
    outfile = f"OT_path_rescaled_DATA{i}"

    curr_sys = ase.io.read(infile, format="lammps-data")
    scaled_pos = curr_sys.get_scaled_positions()
    zvals = curr_sys.numbers

    new_sys = Atoms(scaled_positions = scaled_pos,
                    cell             = relaxed_cell,
                    numbers          = zvals,
                    pbc              = [True, True, True])

    ase.io.write(outfile, new_sys, format="lammps-data", specorder=["Hf","O"], masses=True)
