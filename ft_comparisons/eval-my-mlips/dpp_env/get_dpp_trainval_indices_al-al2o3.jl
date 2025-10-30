using PotentialLearning
using InteratomicPotentials 
using Unitful

xyz_file = "/Users/swyant/cesmix/datasets/al-al2o3_phd/reg_xyz_trainval_no-rd/reg_trainval_no-rd.xyz"
batch_size = 200
ds = load_data(xyz_file, ExtXYZ(u"eV", u"Å")) 


ace = ACE(;species           = [:Al, :O],
          body_order        = 4,
          polynomial_degree = 8,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.15,
          rcutoff           = 5.2)
@show length(ace)

e_descrs = compute_local_descriptors(ds,ace)
confs_kDPP = DataSet(ds .+ e_descrs)
dataset_selector = kDPP(confs_kDPP, GlobalMean(), DotProduct(), batch_size=batch_size)
trainval_inds = get_random_subset(dataset_selector)

open("al-al2o3_dpp_trainval_indices_200.txt", "w") do io
    for idx in trainval_inds
        println(io, idx-1)
    end
end
println("Training indices saved to train_indices.txt")

## Step 2: Create reduced dataset excluding training indices
## Create a mapping from new indices to original indices
#all_inds = collect(1:length(confs_kDPP))
#remaining_inds = setdiff(all_inds, train_inds)  # Indices not in training set
#new_to_original_map = Dict(i => remaining_inds[i] for i in 1:length(remaining_inds))
#
## Create reduced dataset
#reduced_confs = confs_kDPP[remaining_inds]
#
## Step 3: Select validation set from reduced dataset
#val_dataset_selector = kDPP(reduced_confs, GlobalMean(), DotProduct(), batch_size=90)
#val_inds_reduced = get_random_subset(val_dataset_selector)
#
## Step 4: Map validation indices back to original indices
#val_inds_original = [new_to_original_map[i] for i in val_inds_reduced]
#
## Save validation indices to file
#open("new_dpp_validation_indices.txt", "w") do io
#    for idx in val_inds_original
#        println(io, idx-1)
#    end
#end
#println("Validation indices saved to validation_indices.txt")
#
## Summary
#println("\nSummary:")
#println("Total configurations: $(length(confs_kDPP))")
#println("Training set size: $(length(train_inds))")
#println("Remaining configurations: $(length(remaining_inds))")
#println("Validation set size: $(length(val_inds_original))")
