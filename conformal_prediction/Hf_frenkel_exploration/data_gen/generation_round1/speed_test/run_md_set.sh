#!/bin/bash 
mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
