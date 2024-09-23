#!/usr/bin/env python 
from ase.db import connect
from ase import Atoms 
import ase.io
import numpy as np
import pprint

pp = pprint.PrettyPrinter()
large_cell = [[30.,0.,0.],[0.,30.,0.],[0.,0.,30.]]

def convert_row_to_atoms(row):
    atoms = row.toatoms()
    atoms.cell = large_cell
    atoms.pbc = [True,True,True]
    atoms.new_array("forces", np.array(row.data["atomic_forces"],dtype=np.float64))
    atoms.info["energy"] = row.total_energy
    return atoms

def normalized_order(atoms1,atoms2):
    atom1_pos = atoms1.positions
    atom2_pos = atoms2.positions

    #numbers = atoms.numbers
    #annotated_pos = [[numbers[i], *pos[i]] for i in range(len(pos))]
    o_atoms1 = [i for i in range(len(atoms1.numbers)) if atoms1.numbers[i] == 8]
    c_atoms1 = [i for i in range(len(atoms1.numbers)) if atoms1.numbers[i] == 6]
    h_atoms1 = [i for i in range(len(atoms1.numbers)) if atoms1.numbers[i] == 1]

    o_atoms2 = [i for i in range(len(atoms2.numbers)) if atoms2.numbers[i] == 8]
    c_atoms2 = [i for i in range(len(atoms2.numbers)) if atoms2.numbers[i] == 6]
    h_atoms2 = [i for i in range(len(atoms2.numbers)) if atoms2.numbers[i] == 1]
    
    h_rij_list = {}
    for idx1 in h_atoms1:
        h_rij_list[idx1] = []
        for idx2 in h_atoms2:
            rij = np.linalg.norm(atom1_pos[idx1] - atom2_pos[idx2])
            h_rij_list[idx1].append((idx2,rij))
        h_rij_list[idx1] = sorted(h_rij_list[idx1], key=lambda x: x[1])
    
    h_compare_with = [h_rijs[0][0] for k, h_rijs in h_rij_list.items()]

    if len(h_compare_with) != len(set(h_compare_with)):
        jdict = {}
        for atom_key in h_rij_list:
            j = h_rij_list[atom_key][0][0]
            if j in jdict:
                jdict[j].append(atom_key) 
            else:
                jdict[j] = [atom_key]
        #for j in jdict:
        #    assert len(jdict[j]) <3
            
    #pp.pprint(h_rij_list)
    #assert len(h_compare_with) == len(set(h_compare_with))


    o_rij_list = {}
    for idx1 in o_atoms1:
        o_rij_list[idx1] = []
        for idx2 in o_atoms2:
            rij = np.linalg.norm(atom1_pos[idx1] - atom2_pos[idx2])
            o_rij_list[idx1].append((idx2,rij))
        o_rij_list[idx1] = sorted(o_rij_list[idx1], key=lambda x: x[1])

    o_compare_with = [o_rijs[0][0] for k, o_rijs in o_rij_list.items()]

    if len(o_compare_with) != len(set(o_compare_with)):
        jdict = {}
        for atom_key in o_rij_list:
            j = o_rij_list[atom_key][0][0]
            if j in jdict:
                jdict[j].append(atom_key) 
            else:
                jdict[j] = [atom_key]
        #for j in jdict:
        #    assert len(jdict[j]) <3

    #pp.pprint(o_rij_list)

    #assert len(o_compare_with) == len(set(o_compare_with))


    c_rij_list = {}
    for idx1 in c_atoms1:
        c_rij_list[idx1] = []
        for idx2 in c_atoms2:
            rij = np.linalg.norm(atom1_pos[idx1] - atom2_pos[idx2])
            c_rij_list[idx1].append((idx2,rij))
        c_rij_list[idx1] = sorted(c_rij_list[idx1], key=lambda x: x[1])


    c_compare_with = [c_rijs[0][0] for k, c_rijs in c_rij_list.items()]

    if len(c_compare_with) != len(set(c_compare_with)):
        jdict = {}
        for atom_key in c_rij_list:
            j = c_rij_list[atom_key][0][0]
            if j in jdict:
                jdict[j].append(atom_key) 
            else:
                jdict[j] = [atom_key]
        #for j in jdict:
        #    assert len(jdict[j]) <3


    #pp.pprint(c_rij_list)
#   assert len(c_compare_with) == len(set(c_compare_with))


        
    #sorted_indices = [i[0] for i in sorted(enumerate(annotated_pos), key=lambda x:x[1])] # probably won't actually compute consistent order
    #return 
    


# Note that train/test split is between molecules, not MD images

traineq_db = './iso17/reference_eq.db'
eq_images = []
with connect(traineq_db) as conn:
    for row in conn.select():
        atoms_obj = convert_row_to_atoms(row)
        eq_images.append(atoms_obj)

train_images = []
for dbname in ['./iso17/reference.db','./iso17/test_within.db']:
    with connect(dbname) as conn:
        for row in conn.select():
            atoms_obj = convert_row_to_atoms(row)
            train_images.append(atoms_obj)


structure_indices = {}
for i in range(101):
    structure_indices[i] = []
    first_start = i*4000
    #for k in range(4000):
    for k in range(3999):
        idx = first_start + k
        structure_indices[i].append(idx)
    second_start = 404000 + i*1000
    for k in range(1000):
        idx = second_start + k
        structure_indices[i].append(idx)


all_symbols = sorted(list(set([str(image.symbols) for image in train_images]))) 
ase.io.write("my_iso17_trainset.xyz", train_images, format="extxyz")
