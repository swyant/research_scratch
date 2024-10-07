#!/usr/bin/env python 
import ase.io

infile  = "./relaxed_Hf.lammpstrj" 
outfile = "hexag_Hf_round1_DATA"

atoms = ase.io.read(infile, format="lammps-dump-text")
atoms.set_chemical_symbols(["Hf" for _ in range(len(atoms.numbers))])

ase.io.write(outfile, atoms, format="lammps-data", specorder=["Hf", "O"], masses=True)
