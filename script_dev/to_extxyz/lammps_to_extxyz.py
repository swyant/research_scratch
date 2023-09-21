#!/usr/bin/env python
import re 
import sys 
import os.path
import pprint

pp = pprint.PrettyPrinter(indent=2)

def parse_bbox(xlat_toks,ylat_toks,zlat_toks,triclinic):
    if not triclinic:
        assert len(xlat_toks) == len(ylat_toks) == len(zlat_toks) == 2
        (xlo,xhi) = xlat_toks
        (ylo,yhi) = ylat_toks
        (zlo,zhi) = zlat_toks

        latvecs = [ [xhi-xlo,0.0,0.0], [0.0,yhi-ylo,0.0], [0.0,0.0,zhi-zlo] ]
    else:
        assert len(xlat_toks) == len(ylat_toks) == len(zlat_toks) == 3

        # https://docs.lammps.org/Howto_triclinic.html
        (xlob,xhib,xy) = xlat_toks
        (ylob,yhib,xz) = ylat_toks
        (zlob,zhib,yz) = zlat_toks

        xlo = xlob - min(0.0,xy,xz,xy+xz)
        xhi = xhib - max(0.0,xy,xz,xy+xz)
        ylo = ylob - min(0.0,yz)
        yhi = yhib - max(0.0,yz)
        zlo = zlob
        zhi = zhib

        latvecs = [ [xhi-xlo,0.0,0.0], [xy,yhi-ylo,0.0], [xz,yz,zhi-zlo] ]
  
    return latvecs

def parse_dump_forces(dump_fname, all_configs, type_dict):
    with open(dump_fname, "r") as dumpf:
        for line in dumpf:
            if "ITEM: TIMESTEP" in line: 
                config_dict = {}
    
                tstep = int(next(dumpf).strip())
                assert tstep not in all_configs # assume no repeated timesteps, to help ensure consistency between thermo data and dump data

                assert "ITEM: NUMBER OF ATOMS" in next(dumpf)
                natoms = int(next(dumpf).strip())
                config_dict["natoms"] = natoms
    
                # assuming periodic inputs for now
                bbox_header = next(dumpf)
                assert "ITEM: BOX BOUNDS", "pp pp pp:" in bbox_header
    
                if "xy xz yz" in bbox_header:
                    triclinic = True
                else:
                    triclinic = False
    
                xlat_toks = list(map(float,re.split(r"\s+",next(dumpf).strip())))
                ylat_toks = list(map(float,re.split(r"\s+",next(dumpf).strip())))
                zlat_toks = list(map(float,re.split(r"\s+",next(dumpf).strip())))
    
                config_dict["latvecs"] = parse_bbox(xlat_toks,ylat_toks,zlat_toks,triclinic)
    
                # expect very specific information in specific order, could generalize as needed
                assert "ITEM: ATOMS id type x y z fx fy fz" in next(dumpf)
                
                config_species  = []
                config_pos    = []
                config_forces = []
                for _ in range(natoms):
                    atom_toks = re.split(r"\s+",next(dumpf).strip())
                    
                    type_int = int(atom_toks[1])
                    assert type_int in type_map

                    config_species.append(type_map[type_int])
                    config_pos.append(list(map(float,atom_toks[2:5])))
                    config_forces.append(list(map(float,atom_toks[5:8])))
    
                config_dict["species"]  = config_species
                config_dict["pos"]    = config_pos
                config_dict["forces"] = config_forces

                all_configs[tstep] = config_dict

def parse_thermo(thermo_fname, all_configs):
    with open(thermo_fname, "r") as thermof:
        assert "step pe pxx pyy pzz pxy pxz pyz" in next(thermof)
        for line in thermof:
            thermo_toks = re.split(r"\s+",line.strip())
            tstep = int(thermo_toks[0])

            if tstep not in all_configs:
                continue
            else:
                all_configs[tstep]["pe"] = float(thermo_toks[1])

                (pxx,pyy,pzz,pxy,pxz,pyz) = list(map(float,thermo_toks[2:8]))
                all_configs[tstep]["stress"] = [ [pxx,pxy,pxz], [pxy,pyy,pyz], [pxz,pyz,pzz] ]

def write_extxyz(out_fname, all_configs):
    with open(out_fname, "w") as outf:
        for tstep,config in all_configs.items():
            outf.write("{:d}\n".format(config["natoms"]))
            lattice_str = ("Lattice=\"{:.16f} {:.16f} {:.16f} {:.16f} {:.16f} {:.16f} "
                          "{:.16f} {:.16f} {:.16f}\" ").format(*config["latvecs"][0], \
                           *config["latvecs"][1], *config["latvecs"][2])

            property_str = "Properties=species:S:1:pos:R:3:forces:R:3 "

            stress_str = ("stress=\"{:.16f} {:.16f} {:.16f} {:.16f} {:.16f} {:.16f} "
                           "{:.16f} {:.16f} {:.16f}\" ").format(*config["stress"][0], \
                           *config["stress"][1], *config["stress"][2])
            pbc_energy_str = "pbc=\"T T T\" energy={:.16f}\n".format(config["pe"])

            info_line = lattice_str + property_str + stress_str + pbc_energy_str
            outf.write(info_line)

            for i in range(config["natoms"]):
                outf.write("{:} {:32.27f} {:32.27f} {:32.27f} {:32.27f} {:32.27f} {:32.27f}\n".format(\
                            config["species"][i],*config["pos"][i], *config["forces"][i]))




assert len(sys.argv) > 2 # dir type_map1 type_map2...
dirname = str(sys.argv[1])

type_map = {}
for i in range(2,len(sys.argv)):
    type_map[i-1] = sys.argv[i]

#print(type_map)

dump_fname   = os.path.join(dirname, "dump_forces.custom")
thermo_fname = os.path.join(dirname, "thermo.dat") 
assert os.path.exists(dump_fname) and os.path.exists(thermo_fname)

all_configs = {}
parse_dump_forces(dump_fname,all_configs, type_map)
parse_thermo(thermo_fname,all_configs)
#pp.pprint(all_configs)
out_fname = os.path.join(dirname, "configs.xyz")
write_extxyz(out_fname,all_configs)
