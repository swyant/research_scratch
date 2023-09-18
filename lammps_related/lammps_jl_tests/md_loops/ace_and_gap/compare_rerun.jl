include("../../utilities/parse_dump/parse_dump.jl")

configs_orig = fragile_parse_dump2("./dump_gap.custom")
configs_rerun = fragile_parse_dump2("./gap_rerun/dump_rerun_gap.custom")

orig_forces = reshape(reduce(hcat,[config["forces"] for config in configs_orig]), size(configs_orig[1]["forces"])..., :) 
rerun_forces = reshape(reduce(hcat,[config["forces"] for config in configs_rerun]), size(configs_rerun[1]["forces"])..., :) 

@show max_force_discrep = maximum(abs.(rerun_forces -orig_forces)) # 2.0232304720479988e-10
