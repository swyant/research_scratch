#!/bin/zsh

# Pristine runs
for i in {1..10}; do 
  echo "pristine $i"
  python lammps_output2xyz.py "pristine_300K/run_$i" "dump_forces_$i.custom" "pristine_$i.xyz"
done

# Frenkel runs
for i in {1..10}; do 
  echo "Frenkel $i"
  python lammps_output2xyz.py "frenkel_300K/run_$i" "dump_forces_$i.custom" "frenkel_$i.xyz"
done
