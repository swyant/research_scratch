#!/usr/bin/env python
import  glob 
import re
import numpy as np
from ase import Atoms
import ase.io
from copy import deepcopy

H2eV = 27.211386245981
zdict = {'H': 1,
         'C': 6,
         'N': 7,
         'O': 8,
         'F': 9}

large_cell = [[30.,0.,0.],[0.,30.,0.],[0.,0.,30.]]

def parse_atomrefs(fname):
    with open(fname, "r") as f:
        lines = f.readlines()
    ddict = {}
    for line in lines[5:10]:
        tokens = re.split(r"\s+", line.strip())
        elem = tokens[0]
        u0 = H2eV*float(tokens[2])
        u  = H2eV*float(tokens[3])
        h  = H2eV*float(tokens[4])
        g  = H2eV*float(tokens[5])
        ddict[elem] = {"U0": u0,
                       "U" : u,
                       "H" : h,
                       "G" : g}
    return ddict

def parse_uncharacterized(fname):
    with open(fname, "r") as f: 
        lines = f.readlines()
    bad_config_indices = []
    for line in lines[9:3063]:
        index = int(re.split(r"\s+",line.strip())[0])
        bad_config_indices.append(index)
    return bad_config_indices

def modify_bad_float(match_obj):
    lead_coeff = match_obj.group(1)
    remainder  = match_obj.group(2)
    exponent   = match_obj.group(3)
    alt_str    = f"{lead_coeff}.{remainder}e{exponent}"
    return alt_str

atomic_energies = parse_atomrefs("atomref.txt")
bad_indices = parse_uncharacterized("uncharacterized.txt")


xyz_files = glob.glob("./dsgdb9nsd_*.xyz")

file_indices = sorted([re.match(r".*dsgdb9nsd_(\d+)\.xyz$", fname).group(1) for fname in xyz_files])

images = []
for fidx in file_indices:
#fidx = file_indices[66941]
    if int(fidx) in bad_indices:
        continue
    with open(f"./dsgdb9nsd_{fidx}.xyz", "r") as f:
        lines = f.readlines()
    
    num_atoms = int(re.split(r"\s+",lines[0].strip())[0])
    
    line2 = re.split(r"\s+", lines[1].strip())
     
    info_dict = {"gdb_id"            : int(line2[1]),
                 "rotA"              : float(line2[2]),
                 "rotB"              : float(line2[3]),
                 "rotC"              : float(line2[4]),
                 "dipole_moment"     : float(line2[5]),
                 "iso_polarizability": float(line2[6]),
                 "homo"              : H2eV*float(line2[7]),
                 "lumo"              : H2eV*float(line2[8]),
                 "gap"               : H2eV*float(line2[9]),
                 "e_spatial_extent"  : float(line2[10]),
                 "zpve"              : H2eV*float(line2[11]),
                 "U0"                : H2eV*float(line2[12]),
                 "U"                 : H2eV*float(line2[13]),
                 "H"                 : H2eV*float(line2[14]), 
                 "G"                 : H2eV*float(line2[15]),
                 "Cv"                : H2eV*float(line2[16])}
    
    elems = []
    pos   = []
    mull_charges = []
    for i in range(num_atoms):
        atom_toks = re.split(r"\s+",lines[2+i].strip())
        elems.append(atom_toks[0])
        try:
            atom_pos = list(map(float,atom_toks[1:4]))
        except:
            atom_pos_toks = atom_toks[1:4]
            for i in range(3):
                bad_num_match = re.match(r"(-?\d+)\.(\d*)\*\^(-?\d+)$",atom_pos_toks[i])
                if bad_num_match:
                    alt_str = modify_bad_float(bad_num_match)
                    #print(f'bad float for gdb {info_dict["gdb_id"]}, replaced {atom_pos_toks[i]} with {alt_str}')
                    atom_pos_toks[i] = alt_str
            atom_pos = list(map(float, atom_pos_toks))

        pos.append(atom_pos)
        try:
            mull_charge = float(atom_toks[4])
        except:
            bad_num_match = re.match(r"(-?\d+)\.(\d*)\*\^(-?\d+)$",atom_toks[4])
            if bad_num_match:
                mull_charge = modify_bad_float(bad_num_match)
            else:
                print("should error")
        mull_charges.append(mull_charge)
        
    zvals = [zdict[elem] for elem in elems]
    
    # process atomization energies
    atomized_vals = {"U0" : -1*info_dict["U0"],
                     "U"  : -1*info_dict["U"],
                     "H"  : -1*info_dict["H"],
                     "G"  : -1*info_dict["G"]}
    
    for key in ["U0", "U", "H", "G"]:
        for elem in elems:
            atomized_vals[key] += atomic_energies[elem][key]
    
    info_dict["U0_atomized"] = atomized_vals["U0"]
    info_dict["U_atomized"] = atomized_vals["U"]
    info_dict["H_atomized"] = atomized_vals["H"]
    info_dict["G_atomized"] = atomized_vals["G"]
    
    info_dict["energy"] = -1*info_dict["U0_atomized"] # Primary Target
    
    # Raw electronic energy and corresponding formation energy
    info_dict["raw_energy"] = info_dict["U0"] - info_dict["zpve"]
    info_dict["formation_energy"] = -1*info_dict["U0_atomized"] - info_dict["zpve"]
    
    # remaining lines 
    info_dict["vib_freqs"] = np.array(list(map(float,re.split(r"\s+",lines[2+num_atoms].strip()))))
    
    smiles_line = re.split(r"\s+", lines[3+num_atoms].strip())
    info_dict["GDB9_SMILES"] = smiles_line[0]
    info_dict["relaxed_SMILES"] = smiles_line[1]
    
    inchi_line = re.split(r"\s+", lines[4+num_atoms].strip())
    info_dict["GDB9_InChI"]  =  re.match(r"InChI=(.*)$",inchi_line[0]).group(1)
    info_dict["relaxed_InChI"]  = re.match(r"InChI=(.*)$",inchi_line[1]).group(1)
    
    config = Atoms(numbers   = zvals,
                   positions = pos,
                   cell      = deepcopy(large_cell),
                   info      = info_dict,
                   pbc       = [1,1,1])
    
    config.new_array("mulliken_partial_charges", np.array(mull_charges))

    images.append(config)

ase.io.write("qm9_fullset_alldata.xyz", images, format="extxyz")
