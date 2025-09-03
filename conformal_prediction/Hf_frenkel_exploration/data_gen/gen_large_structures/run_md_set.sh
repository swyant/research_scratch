#!/bin/bash 

for i in {3..5}; do
    cd large_8x_frenkel/run_$i 
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../large_dilute_frenkel/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../large_pristine/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../
done

for i in {3..5}; do 
    cd 18x_frenkel/run_$i 
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../18x_dilute_frenkel/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../18x_large_pristine/run_$i
    mpirun -np 20 ~/cesmix/dev/cesmix_lammps/build_fixedPOD_030625/lmp -in nvt.in | tee out.txt
    cd ../../
done
