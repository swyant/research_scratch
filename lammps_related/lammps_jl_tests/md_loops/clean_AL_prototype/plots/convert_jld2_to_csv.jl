using Tables, DataFrames, CSV
using JLD2


all_thermo_data = load("iterative_npt_melt_example_thermo_data.jld2", "all_thermo_data")
all_uq_data = load("iterative_npt_melt_example_uq_data.jld2", "all_uq_data")

for i in 1:length(all_thermo_data)
   num_steps = length(all_thermo_data[i]["all_temps"]) 
   @assert num_steps == length(all_uq_data[i]["energy_stdevs"])

   steps = [idx for idx in 1:num_steps]

   iter_matrix = [steps all_thermo_data[i]["all_temps"] all_thermo_data[i]["all_pes"] all_uq_data[i]["energy_stdevs"]]
   
   CSV.write("iteration$(i)_data_ext.csv", DataFrame(Tables.table(iter_matrix; header=["Steps","Temperature","Potential Energy", "Energy Stdev"])), header=true) 
end
