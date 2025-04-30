#!/bin/bash 

for i in {11..15}; do
    cd frenkel_300K/run_$i 
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../pristine_300K/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../
done
