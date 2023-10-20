#using Distributed; addprocs(4)


include("../../utilities/custom_export2lammps/custom_export2lammps.jl")
include("./train_utils.jl")
include("./run_committee_md_melt.jl")

ds_path = "../../../../../../datasets/HfO2_Sivaraman/prl_2021/raw/train.xyz"
raw_data = read_extxyz(ds_path)
label = "HfO2_N2_P6_r4"

rpib = ACE1.rpi_basis(;
            species = [:Hf, :O],
            N       = 2, 
            maxdeg  = 6,
            D       = ACE1.SparsePSHDegree(; wL = 1.5, csp = 1.0),
            r0      = 2.15,
            rin     = 0.65*2.15,  # Does this meaningfully get used in pin=0?
            rcut    = 4.0,
            pin     = 0,
)

weights = Dict("default" => Dict("E" =>1.0,"F" => 1.0, "V"=>0.0))
vref = JuLIP.OneBody([:Hf => -2.70516846, :O => -0.01277342]...)

initial_data = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in raw_data]



all_new_data = []

label = "HfO2_N2_P6_r4"
yace_files = ["./fits/$(label)_fit$(fit_idx).yace" for fit_idx in 1:5]

yace_lmps = [initialize_committee_member(yace_files[i],2) for i in 2:length(yace_files)] # hard-coding number of types currently...


uq_ddict, thermo_ddict, uncertain_configs  = ace_committee_expts(yace_files[1],yace_lmps; num_steps=300000, vel_seed =12280329, start_temp=3200, end_temp=3200);

new_cfg_tstep = 1
for new_cfg in uncertain_configs
    dummy_cfg = get_gap_eandf(new_cfg["box_bounds"],new_cfg["types"],new_cfg["positions"],new_cfg["masses"],new_cfg_tstep)
    new_cfg_tstep += 1
end

run(`/opt/homebrew/Caskroom/miniforge/base/envs/ase_env/bin/python mod_lammps_to_extxyz.py . Hf O`)

rm("dump_forces.custom", force=true)
rm("thermo.dat", force=true)

for c_lmps in yace_lmps
    LAMMPS.API.lammps_close(c_lmps)
end

cmte_orig_indices = [vec(Matrix(CSV.read("./indices/fit_$(i)_train_idxs.csv",DataFrame,header=false))) for i in 1:5]
new_data = read_extxyz("./new_configs.xyz")

new_data_atoms = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                    weights = weights, v_ref=vref) for at in new_data] 

all_new_data = cat(all_new_data,new_data_atoms,dims=1)
for fit_idx in 1:5

    #train_ds = cat(initial_data[cmte_orig_indices[fit_idx]],all_new_data,dims=1)

    train_ds = initial_data[cmte_orig_indices[fit_idx]]
    #cat didn't preserve eltype? idk
    for adata in all_new_data
        push!(train_ds,adata)
    end

    A, Y, W = ACEfit.assemble(train_ds, rpib)
    solver_reg = ACEfit.QR(; lambda=1e-3)
    results_reg = ACEfit.solve(solver_reg,A,Y)

    pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
    pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
    custom_export2lammps("./new_fits/$(label)_fit$(fit_idx).yace", pot_reg,rpib)
end
#all_train_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
#val_idxs = vec(Matrix(CSV.read("./indices/Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))
#
#perform_and_record_fits(all_data,rpib,vref,label,5;train_subset=all_train_idxs, val_subset=val_idxs)


