#!/bin/zsh

sleep 600 
/Users/swyant/cesmix/dev/cesmix_lammps/build_kokkospod_032124/lmp -in in1.lammps > out1.txt
sleep 5
/Users/swyant/cesmix/dev/cesmix_lammps/build_kokkospod_032124/lmp -in in2.lammps > out2.txt
sleep 5
/Users/swyant/cesmix/dev/cesmix_lammps/build_kokkospod_032124/lmp -in in3.lammps > out3.txt
