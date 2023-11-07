#!/usr/bin/env python
import re 
import numpy as np
import pprint

pp = pprint.PrettyPrinter(indent=1)

base_dump_fname = "dump_unstable_melt.custom"
fcomp_stdevs_fname = "unstable_fcomp_stdevs.csv"

out_dump_fname = "dump_unstable_melt_committee.custom"
fcomp_stdevs = np.loadtxt(fcomp_stdevs_fname)

new_lines = []
with open(base_dump_fname, "r") as dumpf:
    atom_idx = 0 
    for line in dumpf:
        new_lines.append(line)
        if "ITEM: TIMESTEP" in line:  
            new_lines.append(next(dumpf)) 

            lne = next(dumpf)
            assert "ITEM: NUMBER OF ATOMS" in lne
            new_lines.append(lne)

            natoms_lne = next(dumpf) 
            natoms = int(natoms_lne.strip())
            new_lines.append(natoms_lne)
    
            # assuming periodic inputs for now
            bbox_header = next(dumpf)
            new_lines.append(bbox_header)
            assert "ITEM: BOX BOUNDS", "pp pp pp:" in bbox_header
    
            if "xy xz yz" in bbox_header:
                triclinic = True
            else:
                triclinic = False
            
            box_line = next(dumpf) 
            new_lines.append(box_line)
            xlat_toks = list(map(float,re.split(r"\s+",box_line.strip())))
            box_line = next(dumpf) 
            new_lines.append(box_line)
            ylat_toks = list(map(float,re.split(r"\s+",box_line.strip())))
            box_line = next(dumpf) 
            new_lines.append(box_line)
            zlat_toks = list(map(float,re.split(r"\s+",box_line.strip())))
        
            # expect very specific information in specific order, could generalize as needed
            info_line = next(dumpf)
            assert "ITEM: ATOMS id type x y z fx fy fz" == info_line.strip()
            info_line = info_line.strip() + " fx_stdev fy_stdev fz_stdev\n"
            new_lines.append(info_line)
            
            for _ in range(natoms):
                atom_line = next(dumpf)
                atom_line = atom_line.rstrip() + " {:19.14f} {:19.14f} {:19.14f}\n".format(*fcomp_stdevs[atom_idx])
                new_lines.append(atom_line)
                atom_idx += 1 


with open(out_dump_fname, "w") as outdumpf:
    outdumpf.writelines(new_lines)
  
