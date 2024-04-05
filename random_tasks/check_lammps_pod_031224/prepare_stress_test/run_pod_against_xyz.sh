#!/bin/bash 

# hard-coded: number of data files, so manually adjust the for loop and the cat commands

#mkdir -p data_files
#python convert_extxyz_to_lammps.py

mkdir -p ./output
for idx in {0..13}; do 
  echo $idx
  sed -i .bk -e "s|\(log  \).*$|\1\.\/output\/log_$idx|"\
             -e "s|\(read_data     \).*$|\1\.\/data_files\/lammps_pod_stress_test_DATA_$idx|"\
             -e "s|\(variable thermo_fname  string \).*$|\1\.\/output\/thermo_$idx.dat|"\
             -e "s|\(variable dump_fname    string \).*$|\1\.\/output\/dump_forces_$idx.custom|"\
             -e "s|\(variable tstep         equal \).*$|\1$idx|" in.lammps
  ~/cesmix/dev/cesmix_lammps/build_eapod_022024/lmp -in in.lammps >> ./output/out.txt
done

cat output/thermo_0.dat > thermo.dat && awk 'FNR-1' output/thermo_{1..13}.dat >> thermo.dat
cat ./output/dump_forces_{0..13}.custom > ./dump_forces.custom
cat ./output/log_{0..13} > ./logfile

python ../../../script_dev/to_extxyz/lammps_to_extxyz.py ./ Hf O

#zip -r data_files.zip data_files/
#zip -r output.zip output/
