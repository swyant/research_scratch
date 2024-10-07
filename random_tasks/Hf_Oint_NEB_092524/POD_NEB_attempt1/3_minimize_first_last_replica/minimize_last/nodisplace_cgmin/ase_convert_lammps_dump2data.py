#!/usr/bin/env python 
import ase.io

infile  = "./relaxed_HfOint.lammpstrj" 
outfile = "hexag_Hf_tetOint_DATA"

spec_map = {1: "Hf",
            2: "O"}

atoms = ase.io.read(infile, format="lammps-dump-text")
atoms.set_chemical_symbols([spec_map[number] for number in atoms.numbers])
#atoms.set_chemical_symbols(["Hf" for _ in range(len(atoms.numbers))])

ase.io.write(outfile, atoms, format="lammps-data", specorder=["Hf", "O"], masses=True)
