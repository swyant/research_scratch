using PrettyTables


function group_type(d::AtomsData; group_key="config_type")
    group_type = "default"
    for (k,v) in d.atoms.data
        if (lowercase(k)==group_key)
            group_type = v.data
        end
    end
    return group_type
end

function alt_linear_errors(data::AbstractArray{AtomsData}, model; energy_per_atom=true, group_key="config_type", verbose=true)

   mae = Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)
   rmse = Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)
   num = Dict("E"=>0, "F"=>0, "V"=>0)

   config_types = []
   config_mae = Dict{String,Any}()
   config_rmse = Dict{String,Any}()
   config_num = Dict{String,Any}()

   for d in data

       c_t = group_type(d; group_key)
       if !(c_t in config_types)
          push!(config_types, c_t)
          merge!(config_rmse, Dict(c_t=>Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)))
          merge!(config_mae, Dict(c_t=>Dict("E"=>0.0, "F"=>0.0, "V"=>0.0)))
          merge!(config_num, Dict(c_t=>Dict("E"=>0, "F"=>0, "V"=>0)))
       end

       # energy errors
       if !isnothing(d.energy_key)
           if energy_per_atom
               estim = energy(model, d.atoms) / length(d.atoms)
               exact = d.atoms.data[d.energy_key].data / length(d.atoms)
           else
               estim = energy(model, d.atoms)
               exact = d.atoms.data[d.energy_key].data
           end
           mae["E"] += abs(estim-exact)
           rmse["E"] += (estim-exact)^2
           num["E"] += 1
           config_mae[c_t]["E"] += abs(estim-exact)
           config_rmse[c_t]["E"] += (estim-exact)^2
           config_num[c_t]["E"] += 1
       end

       # force errors
       if !isnothing(d.force_key)
           estim = mat(forces(model, d.atoms))
           exact = mat(d.atoms.data[d.force_key].data)
           mae["F"] += sum(abs.(estim-exact))
           rmse["F"] += sum((estim-exact).^2)
           num["F"] += 3*length(d.atoms)
           config_mae[c_t]["F"] += sum(abs.(estim-exact))
           config_rmse[c_t]["F"] += sum((estim-exact).^2)
           config_num[c_t]["F"] += 3*length(d.atoms)
       end

       # virial errors
       if !isnothing(d.virial_key)
           estim = virial(model, d.atoms)[SVector(1,5,9,6,3,2)] ./ length(d.atoms)
           # the following hack is necessary for 3-atom cells:
           #   https://github.com/JuliaMolSim/JuLIP.jl/issues/166
           #exact = d.atoms.data[d.virial_key].data[SVector(1,5,9,6,3,2)] ./ length(d.atoms)
           exact = hcat(d.atoms.data[d.virial_key].data...)[SVector(1,5,9,6,3,2)] ./ length(d.atoms)
           mae["V"] += sum(abs.(estim-exact))
           rmse["V"] += sum((estim-exact).^2)
           num["V"] += 6
           config_mae[c_t]["V"] += sum(abs.(estim-exact))
           config_rmse[c_t]["V"] += sum((estim-exact).^2)
           config_num[c_t]["V"] += 6
       end
    end

    # finalize errors
    for (k,n) in num
        (n==0) && continue
        rmse[k] = sqrt(rmse[k] / n)
        mae[k] = mae[k] / n
    end
    errors = Dict("mae"=>mae, "rmse"=>rmse)

    # finalize config errors
    for c_t in config_types
        for (k,c_n) in config_num[c_t]
            (c_n==0) && continue
            config_rmse[c_t][k] = sqrt(config_rmse[c_t][k] / c_n)
            config_mae[c_t][k] = config_mae[c_t][k] / c_n
        end
    end
    config_errors = Dict("mae"=>config_mae, "rmse"=>config_rmse)

    # merge errors into config_errors and return
    push!(config_types, "set")
    merge!(config_errors["mae"], Dict("set"=>mae))
    merge!(config_errors["rmse"], Dict("set"=>rmse))

    if verbose
        print_errors_tables(config_errors)
    end 

    return config_errors
end


function print_errors_tables(config_errors::Dict)
    print_rmse_table(config_errors)
    print_mae_table(config_errors)
end

function _print_err_tbl(D::Dict)
    header = ["Type", "E [meV]", "F [eV/A]", "V [meV]"]
    config_types = setdiff(collect(keys(D)), ["set",])
    push!(config_types, "set")
    table = hcat(
        config_types,
        [1000*D[c_t]["E"] for c_t in config_types],
        [D[c_t]["F"] for c_t in config_types],
        [1000*D[c_t]["V"] for c_t in config_types],
    )
    pretty_table(
        table; header=header,
        body_hlines=[length(config_types)-1],
        formatters=ft_printf("%5.3f"),
        crop = :horizontal)

end

function print_rmse_table(config_errors::Dict; header=true)
    if header; (@info "RMSE Table"); end 
    _print_err_tbl(config_errors["rmse"])
end

function print_mae_table(config_errors::Dict; header=true)
    if header; (@info "MAE Table"); end 
    _print_err_tbl(config_errors["mae"])
end

