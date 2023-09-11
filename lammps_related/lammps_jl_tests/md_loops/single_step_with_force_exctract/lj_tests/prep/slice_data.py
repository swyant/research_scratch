#!/usr/bin/env python
import re 

input_data_fname = "fcc_lj_Ar_box-zeroed_shifted_DATA"
output_data_fname = "fcc_lj_Ar_smaller_DATA"

zcut = 19.8
zbuff = 1.3186998295724957

parse_atoms = False
orig_latspecs = []
orig_atoms = []
num_types = -1
mass_dict = {}
with open(input_data_fname, "r") as inputf: 
  for line in inputf:
    if parse_atoms:
      atom_tokens = re.split(r"\s+", line.strip())

      atomid   = int(atom_tokens[0])
      atomtype = int(atom_tokens[1])
      xcoord   = float(atom_tokens[2])
      ycoord   = float(atom_tokens[3])
      zcoord   = float(atom_tokens[4])

      orig_atoms.append([atomid,atomtype,xcoord,ycoord,zcoord])

    else:
      atom_toks = re.split(r"\s+", line.strip())
      if atom_toks[-1] == "types":
        num_types = int(atom_toks[0])
        print(num_types)

      if atom_toks[-1] == "xhi":
        latspec = [float(atom_toks[0]), float(atom_toks[1])]
        orig_latspecs.append(latspec)
        
        for _ in range(2):
          atom_toks = re.split(r"\s+", next(inputf).strip())
          latspec = [float(atom_toks[0]), float(atom_toks[1])]
          orig_latspecs.append(latspec)
        print(orig_latspecs)

      if atom_toks[0] == "Masses":
        next(inputf)
        
        mass_line = next(inputf)
        mass_line_toks = re.split(r"\s+", mass_line.strip())
        while not mass_line_toks[0] == "":
          atom_type_num =  int(mass_line_toks[0])
          mass = float(mass_line_toks[1])
          mass_dict[atom_type_num] = mass 
          mass_line = next(inputf)
          mass_line_toks = re.split(r"\s+", mass_line.strip())
        
        print(mass_dict)

      if line.strip() == "Atoms":
        parse_atoms = True
        next(inputf)

new_coords = []
new_id = 1
largest_z = -10000
for orig_atom in orig_atoms:
  if orig_atom[-1] < zcut:
    new_coords.append([new_id, *orig_atom[1:]])
    new_id += 1
    
    if orig_atom[-1] > largest_z:
      largest_z = orig_atom[-1]
#print(new_coords)

num_atoms = len(new_coords)
new_z = largest_z + zbuff

with open(output_data_fname, "w") as outf:
  outf.write("LAMMPS DATA file\n\n")
  outf.write("{:d} atoms\n\n".format(num_atoms))
  outf.write("{:d} atom types\n\n".format(num_types))

  outf.write("{:.16e} {:.16e} xlo xhi\n".format(*orig_latspecs[0]))
  outf.write("{:.16e} {:.16e} ylo yhi\n".format(*orig_latspecs[1]))
  outf.write("{:.16e} {:.16e} zlo zhi\n\n".format(orig_latspecs[2][0], new_z))

  outf.write("Masses\n\n")
  for k,v in mass_dict.items():
    outf.write("{:d} {:f}\n".format(k,v))

  outf.write("\nAtoms\n\n")
  for new_coord in new_coords:
    outf.write("{:>4d} {:1d} {:22.16f} {:22.16f} {:22.16f}\n".format(*new_coord))

