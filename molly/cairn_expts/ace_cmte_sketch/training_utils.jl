using Random

function obtain_train_idxs(frac::Float64, set_size::Int64)
    @assert frac <= 1.0
    num_select = Int(floor(frac*set_size))

    perm_idxs = Random.randperm(set_size)
    rand_set_idxs = perm_idxs[begin:1:num_select]

    rand_set_idxs
end

function obtain_train_idxs(frac::Float64, train_subset::Vector{Int64})
    @assert frac <= 1.0
    
    set_size = length(train_subset)
    num_select = Int(floor(frac*set_size))

    perm_idxs = Random.randperm(set_size)
    rand_set_idxs = perm_idxs[begin:1:num_select]

    train_subset[rand_set_idxs]
end


function perform_and_record_simple_fits(ds, basis, num_fits; frac=0.7, weights=[100.0,1.0], int=false, train_subset=nothing)
    lbps = Vector{LinearBasisPotential}()
    for fit_idx in 1:num_fits
        print("###### RUNNING FIT $(fit_idx) ######")
        if isnothing(train_subset)
            train_idxs = obtain_train_idxs(frac, length(ds))
        else
            train_idxs = obtain_train_idxs(frac, train_subset)
        end
        CSV.write("./initial_fits/fit$(fit_idx)_train_idxs.csv",DataFrame(Tables.table(train_idxs)),header=false)

        lbp = LBasisPotential(basis)
        learn!(lbp,ds[train_idxs],weights,int) 

        CSV.write("./initial_fits/fit$(fit_idx)_coeffs.csv", DataFrame(Tables.table(lbp.β)),header=false)
        push!(lbps,lbp)    
    end

    lbps
end


function perform_and_record_pce_fits(trainset::Vector{<:AbstractSystem}, pce_template::PolynomialChaos, ref, num_fits; frac=0.7, weights=[10000.0,1.0], intcpt=false, train_subset=nothing)
    pces = Vector{PolynomialChaos}()
    for fit_idx in 1:num_fits
        print("RUNNING PCE FIT $(fit_idx)\n")
        if isnothing(train_subset)
            train_idxs = obtain_train_idxs(frac, length(trainset))
        else
            train_idxs = obtain_train_idxs(frac, train_subset)
        end
        CSV.write("./pce_fits/lq_fit$(fit_idx)_train_idxs.csv",DataFrame(Tables.table(train_idxs)),header=false)

        lp = learn!(trainset[train_idxs],ref, pce_template, weights,intcpt) 

        CSV.write("./pce_fits/lq_fit$(fit_idx)_coeffs.csv", DataFrame(Tables.table(lp.β)),header=false)
        pce_entry = deepcopy(pce_template)
        pce_entry.params = lp.β

        push!(pces,pce_entry)    
    end
    pces
end

