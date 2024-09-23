#!/usr/bin/env python 
from ase.db import connect
from ase import Atoms 
import ase.io
import numpy as np

large_cell = [[30.,0.,0.],[0.,30.,0.],[0.,0.,30.]]

def convert_row_to_atoms(row):
    atoms = row.toatoms()
    atoms.cell = large_cell
    atoms.pbc = [True,True,True]
    atoms.new_array("forces", np.array(row.data["atomic_forces"],dtype=np.float64))
    atoms.info["energy"] = row.total_energy
    return atoms

# Note that my train/test split is only between molecules, not MD images

#traineq_db = './iso17/reference_eq.db'
#eq_images = []
#with connect(traineq_db) as conn:
#    for row in conn.select():
#        atoms_obj = convert_row_to_atoms(row)
#        eq_images.append(atoms_obj)

#train_images = []
#for dbname in ['./iso17/reference.db','./iso17/test_within.db']:
#    with connect(dbname) as conn:
#        for row in conn.select():
#            atoms_obj = convert_row_to_atoms(row)
#            train_images.append(atoms_obj)

test_images = []
for dbname in ['iso17/test_other.db']:
    with connect(dbname) as conn:
        for row in conn.select():
            atoms_obj = convert_row_to_atoms(row)
            test_images.append(atoms_obj)


#ase.io.write("my_iso17_train.xyz", train_images, format="extxyz")
ase.io.write("my_iso17_test.xyz", test_images, format="extxyz")
