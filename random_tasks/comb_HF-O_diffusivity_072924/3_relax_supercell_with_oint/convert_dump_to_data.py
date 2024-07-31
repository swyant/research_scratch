#!/usr/bin/env python

import ase.io
import sys

#input_file = "intermediate_Hf_supercell_Oint_nonvt.lammpstrj"
#output_file = "intermediate_Hf_supercell_Oint_nonvt_DATA"

#input_file = "relaxed_Hf_supercell_Oint_nonvt.lammpstrj"
#output_file = "relaxed_orthog_Hf_15_9_9_singleOint_DATA"

input_file = "better_Oint_starting_location/alt_relaxed_Hf_supercell_Oint_nonvt.lammpstrj"
output_file = "better_Oint_starting_location/alt_relaxed_Hf_supercell_Oint_nonvt_DATA"

atoms = ase.io.read(input_file, specorder=["Hf","O"])

ase.io.write(output_file, atoms, format="lammps-data", specorder=["Hf", "O"], masses=True, atom_style="charge")
