using LAMMPS

lmp = LMP(["-screen","none"]) 
command(lmp, "log none")
command(lmp, "read_data fake")
command(lmp, "units          metal")
command(lmp, "atom_style     atomic")
command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step

command(lmp, "read_data     tetrag_hfo2_sample_DATA")

command(lmp, "pair_style     pace")
command(lmp, "pair_coeff     * * dummy.yace Hf O")
