#!/bin/bash 

for n in {1..10};
do 
  ~/local_software/lammps/lammps_081723/clang_shared_build/lmp -in in.lammps >> out.txt 
done
