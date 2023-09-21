#!/usr/bin/env python
import re 
import sys 
import os.path

dirname = str(sys.argv[1])
#print(dirname)

dump_fname   = os.path.join(dirname, "dump_forces.custom")
thermo_fname = os.path.join(dirname, "thermo.dat") 
assert os.path.exists(dump_fname)
assert os.path.exists(thermo_fname)


configs = []
with open(dump_fname, "r") as dumpf:
  for line in dumpf:
    print(line)
    if "ITEM: TIMESTEP" in line: 
        print(next(dumpf))
        print(next(dumpf))
