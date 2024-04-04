#!/usr/bin/env python

import ase.io
import sys 

input_fname = str(sys.argv[1])
output_fname = str(sys.argv[2])

atoms = ase.io.read(input_fname)
ase.io.write(output_fname,atoms,format="lammps-data", masses=True)
