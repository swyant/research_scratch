from ase.io import read as ase_read
from ase.io import write as ase_write
import sys 

input_file = str(sys.argv[1])
output_file = str(sys.argv[2])

atoms = ase_read(input_file, format="lammps-dump-text")

# hard-coded for systems that are just Hf
atoms.set_chemical_symbols(["Hf"]*len(atoms))

ase_write(output_file, atoms, format="extxyz")
