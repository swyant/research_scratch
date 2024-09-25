#!/usr/bin/env python

import numpy as np 
from ase import Atoms
import ase.io
from copy import deepcopy

kcalmol2eVatom = 4184*6.241509074*10**18/(6.02214076*10**23)
large_cell = [[30.,0.,0.],[0.,30.,0.],[0.,0.,30.]]
pbc = [True, True, True]

fnames = {'uracil'       : './rmd17/npz_data/rmd17_uracil.npz',
          'salicylic'    : './rmd17/npz_data/rmd17_salicylic.npz',
          'paracetamol'  : './rmd17/npz_data/rmd17_paracetamol.npz',
          'malonaldehyde': './rmd17/npz_data/rmd17_malonaldehyde.npz',
          'naphthalene'  : './rmd17/npz_data/rmd17_naphthalene.npz',
          'benzene'      : './rmd17/npz_data/rmd17_benzene.npz',
          'aspirin'      : './rmd17/npz_data/rmd17_aspirin.npz',
          'toluene'      : './rmd17/npz_data/rmd17_toluene.npz',
          'azobenzene'   : './rmd17/npz_data/rmd17_azobenzene.npz',
          'ethanol'      : './rmd17/npz_data/rmd17_ethanol.npz'}


for mol in fnames:
    print(f"Molecule: {mol}")

    configs = []
    data = np.load(fnames[mol])

    num_images = data['energies'].size
    zvals = data['nuclear_charges']

    for i in range(num_images):
        if i%100 == 0:
            print(i)

        pos = data['coords'][i,:,:]
        energy = data["energies"][i]*kcalmol2eVatom
        forces = data["forces"][i,:,:]*kcalmol2eVatom
        
        info_dict = {"energy": energy}

        config = Atoms(numbers   = zvals, 
                       positions = pos,
                       cell      = deepcopy(large_cell),
                       info      = info_dict,
                       pbc       = pbc)

        config.new_array("forces", forces)
        configs.append(config)

#    ase.io.write(f"rmd17_{mol}.xyz", configs, format="extxyz")
