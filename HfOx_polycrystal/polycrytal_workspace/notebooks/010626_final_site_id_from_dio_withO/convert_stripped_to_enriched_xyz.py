#!/usr/bin/env python

from ase.io import read as ase_read
from ase.io import write as ase_write
import numpy as np
import copy
import json

stripped_atoms = ase_read("test_DATA", format="lammps-data")

# VERY HACK DOUBLE MAPS BECAUSE THIS IS QUICK AND DIRTY
with open("noOidx2orig.json", "r") as f:
    index_map = json.load(f)
# I want to reverse this, i.e. go from orig to noO
index_map = {int(v): int(k) for k,v in index_map.items()}

with open("test_DATA.hf_mapping.json", "r") as f:
    back2orig_map = json.load(f)
back2orig_map = {int(k): int(v) for k,v, in back2orig_map.items()} # this is in the correct order

grain_ptm_data = np.load("grains_ptm_111025_min4_fixed.npz") # notice! using the fixed version now
noO_grains = grain_ptm_data["grains"]
noO_ptm_types = grain_ptm_data["ptm_types"]

xyz_grain_idxs = []
xyz_ptm_types = []

for i,atm in enumerate(stripped_atoms):
    if atm.symbol == "O":
        xyz_grain_idxs.append(-1)
        xyz_ptm_types.append(-1)
        continue

    xyz_grain_idxs.append(noO_grains[index_map[back2orig_map[i]]])
    xyz_ptm_types.append(noO_ptm_types[index_map[back2orig_map[i]]])

stripped_atoms_out = copy.deepcopy(stripped_atoms)
stripped_atoms_out.set_array("grain_index", np.array(xyz_grain_idxs))
stripped_atoms_out.set_array("ptm_type", np.array(xyz_ptm_types))

ase_write("test_DATA.xyz", stripped_atoms_out, format="extxyz")
