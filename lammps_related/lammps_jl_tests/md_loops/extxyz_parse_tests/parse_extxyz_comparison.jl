using PotentialLearning
using Unitful, UnitfulAtomic 

include("../../utilities/parse_dump/parse_dump.jl")
include("../../utilities/fixed_parse_extXYZ/parse_extXYZ.jl")

configs_man1 = fragile_parse_dump2("./dump_forces1.custom")
configs_man1 = parse_thermo_file("./thermo1.dat"; config_list=configs_man1)

man_force_comps1 = [force_comp for config in configs_man1 for force_comp in config["forces"]]
man_energies1 = [config["pe"] for config in configs_man1]

ds1 = fixed_load_data("./configs1.xyz",ExtXYZ(u"eV",u"Å"))
force_comps1 = get_all_forces(ds1)
energies1 = get_all_energies(ds1)

@show max_force_discrep1 = maximum(abs.(man_force_comps1 - force_comps1)) #0.0
@show max_pe_discrep1 = maximum(abs.(man_energies1-energies1)) #0.0 


### comparing concatenated
configs_man2 = fragile_parse_dump2("./dump_forces2.custom")
configs_man2 = parse_thermo_file("./thermo2.dat"; config_list=configs_man2)

man_force_comps2 = [force_comp for config in configs_man2 for force_comp in config["forces"]]
man_energies2 = [config["pe"] for config in configs_man2]

all_force_comps_man = cat(man_force_comps1,man_force_comps2; dims=1)
all_energies_man = cat(man_energies1,man_energies2; dims=1)

ds2=fixed_load_data("./configs_all.xyz", ExtXYZ(u"eV",u"Å"))

force_comps_all = get_all_forces(ds2)
energies_all = get_all_energies(ds2)

@show max_force_discrep_all = maximum(abs.(all_force_comps_man - force_comps_all)) #0.0
@show max_pe_discrep_all = maximum(abs.(all_energies_man-energies_all)) #0.0


configs_man_nonorthog = fragile_parse_dump2("dump_forces_nonorthog.custom")
configs_man_nonorthog = parse_thermo_file("thermo_nonorthog.dat"; config_list=configs_man_nonorthog)

man_force_comps_nonorthog = [force_comp for config in configs_man_nonorthog for force_comp in config["forces"]]
man_energies_nonorthog = [config["pe"] for config in configs_man_nonorthog]

ds_nonorthog = fixed_load_data("configs_nonorthog.xyz", ExtXYZ(u"eV",u"Å"))

force_comps_nonorthog = get_all_forces(ds_nonorthog)
energies_nonorthog = get_all_energies(ds_nonorthog)

@show max_force_discrep_nonorthog = maximum(abs.(man_force_comps_nonorthog - force_comps_nonorthog)) #0.0
@show max_pe_discrep_all = maximum(abs.(man_energies_nonorthog - energies_nonorthog)) #0.0