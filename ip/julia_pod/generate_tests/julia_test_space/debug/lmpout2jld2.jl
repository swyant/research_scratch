using DelimitedFiles
using JLD2 

dir = "./3body_monoclinic_hfo2"
out = "monoclinic_hfo2_3body_output.jld2"

pe_file = "$(dir)/rerun_pe.dat"
f_file  = "$(dir)/f.txt"

f = readdlm(f_file)
pe = readdlm(pe_file, skipstart=1)[2]

save(out, Dict("pe"  => pe,
               "f"   => f))
