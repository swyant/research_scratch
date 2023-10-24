
include("./single_point_gap.jl")
include("./custom_export2lammps.jl")
#= 
Why do it this way: I may want to naively take a random subset of the full data, or 
 I already have a list of subset idx, and I'm taking an additional subset of that, but 
 I still want the idx to be wrt the original dataset
 =#
function obtain_train_idxs(frac; train_subset=none, set_size=nothing)
    @assert frac <= 1.0
    if !isnothing(train_subset)
        if !isnothing(set_size)
            print("set size will be overwritten")
        end
        set_size = length(train_subset)
    end
    num_select = Int(floor(frac*set_size))

    perm_idxs = Random.randperm(set_size)
    rand_set_idxs = perm_idxs[begin:1:num_select]

    if isnothing(train_subset)
        return rand_set_idxs
    else 
        return train_subset[rand_set_idxs]
    end
end

function perform_and_record_fits(full_ds, rpib, vref, label, num_fits; frac=0.9, train_subset=nothing, val_subset=nothing)
    for fit_idx in 1:num_fits
        print("###### RUNNING FIT $(fit_idx) ######")
        if isnothing(train_subset)
            train_idxs = obtain_train_idxs(frac; set_size=length(full_ds))
        else
            train_idxs = obtain_train_idxs(frac; train_subset=train_subset)
        end

        CSV.write("indices/fit_$(fit_idx)_train_idxs.csv",DataFrame(Tables.table(train_idxs)),header=false)
    
        A, Y, W = ACEfit.assemble(full_ds[train_idxs], rpib)
        solver_reg = ACEfit.QR(; lambda=1e-3)
        results_reg = ACEfit.solve(solver_reg,A,Y)
        CSV.write("fits/$(label)_fit$(fit_idx)_coeffs.csv", DataFrame(Tables.table(results_reg["C"])),header=false)
    
        pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
        pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
        train_errors = ACE1pack.linear_errors(full_ds[train_idxs],pot_reg)
        if !isnothing(val_subset)
            val_errors = ACE1pack.linear_errors(full_ds[val_subset],pot_reg)
        end
    
        save("./errors/$(label)_fit$(fit_idx)_errors.jld2","train_errors",train_errors,"val_errors",val_errors)
        custom_export2lammps("./fits/$(label)_fit$(fit_idx).yace", pot_reg,rpib)
    end
end

function label_new_configs(uq_configs, weights, vref; start_tstep=1)
    tstep = start_tstep
    for new_cfg in uq_configs
        dummy_cfg = get_gap_eandf(new_cfg["box_bounds"],new_cfg["types"],new_cfg["positions"],new_cfg["masses"],tstep)
        tstep += 1
    end

    run(`/opt/homebrew/Caskroom/miniforge/base/envs/ase_env/bin/python new_code/mod_lammps_to_extxyz.py . Hf O`)
    rm("dump_forces.custom", force=true)
    rm("thermo.dat", force=true)

    new_data = read_extxyz("./new_configs.xyz")

    new_data_atoms = [ AtomsData(at;  energy_key = "energy", force_key="forces",
                        weights = weights, v_ref=vref) for at in new_data] 
    
    new_data_atoms, tstep
end 


function updated_ace_committee_fits(initial_data,new_configs,rpib,vref,num_members; initial_indices=nothing)
    @assert length(initial_indices) == num_members
    for fit_idx in 1:num_members   
        train_ds = initial_data[initial_indices[fit_idx]]
        #cat didn't preserve eltype? idk
        for adata in new_configs
            push!(train_ds,adata)
        end
    
        A, Y, W = ACEfit.assemble(train_ds, rpib)
        print("finished assemble")
        solver_reg = ACEfit.QR(; lambda=1e-3)
        results_reg = ACEfit.solve(solver_reg,A,Y)
    
        pot_reg_tmp = JuLIP.MLIPs.combine(rpib,results_reg["C"])
        pot_reg = JuLIP.MLIPs.SumIP(pot_reg_tmp,vref)
        custom_export2lammps("./new_fits/$(label)_fit$(fit_idx).yace", pot_reg,rpib)
    end
end
