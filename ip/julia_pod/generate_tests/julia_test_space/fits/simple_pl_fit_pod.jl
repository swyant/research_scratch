using PotentialLearning, InteratomicPotentials 
using Unitful
using DelimitedFiles

train_file = "../../files/hfo2_training/hfo2_training_set.xyz"
param_file = "../../fits/hfo2_4-body/sample_4body_hfo2_param.pod"
elems = [:Hf, :O]
out_file = "hfo2_lbp_beta.txt"

ds = load_data(train_file, ExtXYZ(u"eV", u"Å"))

lmp_pod = LAMMPS_POD(param_file, elems)
lbp = LBasisPotential(lmp_pod) 

ws, int = [100.0^2, 2.0^2], false
@time _AtWA, _AtWb  = PotentialLearning.ooc_learn!(lbp, ds; 
                                                   ws=ws,
                                                   reg_style=:scale_thresh,
                                                   symmetrize=true, 
                                                   λ=1e-10)

open(out_file, "w") do io
    writedlm(io, lbp.β)
end
