#!/usr/bin/env python
import os 
from copy import deepcopy
import numpy as np 
from apsw import Connection
from ase import Atoms
import ase.io

BOHR2ANG = 0.529177210544 # 2022 CODATA, multiply to go from Bohr 2 angstrom (opposite of what they had in their script)
species_map = {1: 'H',
               6: 'C',
               7: 'N',
               8: 'O',
               9: 'F'}

large_cell = [[30.,0.,0.],[0.,30.,0.],[0.,0.,30.]]

repo_dir = '/home/swyant/cesmix/dev/other/AIRS/OpenDFT/QHBench/QH9'
connection = Connection(os.path.join(repo_dir, 'datasets/QH9Dynamic_300k/QH9_Dyn_300k.db'))
cursor = connection.cursor()

data = cursor.execute("select * from data")
images = []
for row in data:
    pos = np.frombuffer(row[4],np.float64)
    num_atoms = row[2]
    pos = BOHR2ANG*pos.reshape(num_atoms,3)
    
    zvals = np.frombuffer(row[3],np.int32)
    elems = [species_map[z] for z in zvals]
    
    mol_id = int(np.frombuffer(row[0], dtype=np.int64)[0])
    geo_id = row[1]
    tstep  = row[8]
    
    energy = row[6]
    ekin   = row[5]

    info_dict = {"mol_id": mol_id, 
                 "geo_id": geo_id,
                 "tstep" : tstep,
                 "energy": energy, 
                 "ekin"  : ekin}

    config = Atoms(numbers   = zvals,
                   positions = pos,
                   cell      = deepcopy(large_cell),
                   info      = info_dict)
    
    images.append(config)

ase.io.write("full_set.xyz", images, format="extxyz")
