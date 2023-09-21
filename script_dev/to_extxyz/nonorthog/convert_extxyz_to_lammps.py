import ase.io

atoms = ase.io.read("./dump.xyz")
ase.io.write("test_nonorthog_DATA",atoms,format="lammps-data")
