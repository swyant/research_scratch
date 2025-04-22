#!/bin/zsh

# Large 8x Frenkel runs
for i in {1..2}; do 
  echo "large_8x_frenkel $i"
  python lammps_output2xyz.py "large_8x_frenkel/run_$i" "dump_forces_1.custom" "large_8x_frenkel_$i.xyz"
done

# Large Dilute Frenkel runs
for i in {1..2}; do 
  echo "large_dilute_frenkel $i"
  python lammps_output2xyz.py "large_dilute_frenkel/run_$i" "dump_forces_1.custom" "large_dilute_frenkel_$i.xyz"
done

# Large Pristine
for i in {1..2}; do 
  echo "large_pristine $i"
  python lammps_output2xyz.py "large_pristine/run_$i" "dump_forces_1.custom" "large_pristine_$i.xyz"
done

