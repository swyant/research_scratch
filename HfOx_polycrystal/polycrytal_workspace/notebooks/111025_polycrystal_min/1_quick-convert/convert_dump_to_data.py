from ase.io import read as ase_read 
from ase.io import write as ase_write 
import sys 

input_file = str(sys.argv[1])
output_file = str(sys.argv[2]) 

atoms = ase_read(input_file, format="lammps-dump-text")
ase_write(output_file, atoms, format="lammps-data", masses=True, specorder=["Hf", "O"])
