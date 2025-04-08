#!/bin/bash 

for i in {1..2}; do
    cd large_8x_frenkel/run_$i 
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../large_dilute_frenkel/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../large_pristine/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../
done
