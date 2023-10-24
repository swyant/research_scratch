
# load things 
include("./new_code/dependencies.jl")
include("./new_code/md_active_learning.jl")
include("./new_code/train_utils.jl")

# initialization function 
function simple_ace_nvt(;driver_yace_fname="sample.yace",temp=300, vel_seed=12280329, atom_type_list=["Hf","O"]) 
    atomtype_str = ""
    for elem_str in atom_type_list
        atomtype_str = atomtype_str * " " * elem_str
    end


    lmp = LMP(["-screen","none"])
    command(lmp, "log none")

    command(lmp, "units          metal")
    command(lmp, "boundary       p p p")
    command(lmp, "atom_style     atomic")
    command(lmp, "neigh_modify   delay 0 every 1 check no") # neighborlist rebuilt every step
    
    command(lmp, "read_data     tetrag_hfo2_sample_DATA")
    
    command(lmp, "pair_style     pace")
    command(lmp, "pair_coeff     * * $(driver_yace_fname)$(atomtype_str)")
    
    command(lmp, "variable T        equal  $(temp)")
    command(lmp, "variable Tdamp    equal  0.1")
    command(lmp, "variable Tseed    equal  $(vel_seed)")
    #command(lmp, "variable dumpf    equal  100")
    
    command(lmp, "velocity     all create \${T} \${T} mom yes rot yes dist gaussian")
    command(lmp, "fix          nvt all nvt temp \${T} \${T} \${Tdamp}")

    println("successfully ran the initialization commands")

    lmp
end 



# initial parameters 
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

all_train_idxs = vec(Matrix(CSV.read("./Siviraman_HfO2_my_train_idxs.csv",DataFrame,header=false)))
val_idxs = vec(Matrix(CSV.read("./Siviraman_HfO2_my_val_idxs.csv",DataFrame,header=false)))

#initial_committee_fits(all_data,rpib,vref,label,5;train_subset=all_train_idxs, val_subset=val_idxs)
initial_indices = [vec(Matrix(CSV.read("./initial_indices/fit_$(i)_train_idxs.csv",DataFrame,header=false))) for i in 1:5]
smaller_initial_indices= [init_indices[1:500] for init_indices in initial_indices]

all_uq_data = []
all_thermo_data = []
new_configs = []
new_cfg_tstep = 1
for iter_index in 1:5
#iter_index = 1
    if iter_index == 1
        fit_dir = "./initial_fits"
    else
        fit_dir = "./new_fits"
    end
    println(fit_dir)
    
    yace_files = ["$(fit_dir)/$(label)_fit$(fit_idx).yace" for fit_idx in 1:5]
    ace_cmte = PACE_Ensemble(yace_files, ["Hf","O"]; energy_thresh=0.1)
    
    nvt_params = Dict(:driver_yace_fname => yace_files[1],
                      :temp              => 4500,
                      :vel_seed          =>12280329,
                      :atom_type_list    => ["Hf","O"],
                      )

    uq_data, thermo_data, uq_configs = md_activelearn(simple_ace_nvt,nvt_params,ace_cmte,10000) 
    #close!(ace_cmte) # I think this fucks with the garbage collector :(

    push!(all_uq_data,deepcopy(uq_data))
    push!(all_thermo_data, deepcopy(thermo_data))

    if length(uq_configs) > 1
        new_atoms_data, new_cfg_tstep = label_new_configs(uq_configs,weights,vref; start_tstep=new_cfg_tstep)
        for new_atoms in new_atoms_data
            push!(new_configs,new_atoms)
        end

        updated_ace_committee_fits(initial_data,new_configs, rpib,vref,5; initial_indices=initial_indices)
    else
        println("Completed")
        break
    end
end


