# Run multiple POD fitting experiments in serial or parallel using LAMMPS

using IterTools, OrderedCollections


# Parallel execution. Warning: a high number of parallel experiments may degrade system performance.
parallel = false

# path
experiments_path = "./experiments/"
lammps_path = "../../../lammps/build/lmp"

# param.pod ####################################################################

params = OrderedDict()

n = 15

# chemical element symbols
params[:species]                                 = ["Hf O" for _ in 1:n]

# periodic boundary conditions
params[:pbc]                                     = ["1 1 1" for _ in 1:n]

# inner cut-off radius
params[:rin]                                     = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# outer cut-off radius
params[:rcut]                                    = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]

# polynomial degrees for radial basis functions
params[:bessel_polynomial_degree]                = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
params[:inverse_polynomial_degree]               = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

# one-body potential
params[:onebody]                                 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# two-body linear POD potential
params[:twobody_number_radial_basis_functions]   = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]

# three-body linear POD potential
params[:threebody_number_radial_basis_functions] = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
params[:threebody_angular_degree]                = [2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4]

# four-body linear POD potential
params[:fourbody_number_radial_basis_functions]  = [0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
params[:fourbody_angular_degree]                 = [0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4]

params[:true4BodyDesc]                           = [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

# five-body linear POD potential
params[:fivebody_number_radial_basis_functions]  = [0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4]
params[:fivebody_angular_degree]                 = [0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4, 2, 3, 4]

# six-body linear POD potential
params[:sixbody_number_radial_basis_functions]   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3]
params[:sixbody_angular_degree]                  = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 2, 3, 4]

# seven-body linear POD potential
params[:sevenbody_number_radial_basis_functions] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2]
params[:sevenbody_angular_degree]                = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4]

# data.pod #####################################################################

data = OrderedDict()

data[:file_format] = ["extxyz"]
data[:file_extension] = ["extxyz"]

data[:path_to_training_data_set] = ["\"../../train\""]
data[:path_to_test_data_set] = ["\"../../test\""]

data[:fitting_weight_energy] = [100.0]
data[:fitting_weight_force] = [1.0]
data[:fitting_regularization_parameter] = [1e-8]

data[:error_analysis_for_training_data_set] = [1]
data[:error_analysis_for_test_data_set] = [1]

data[:fraction_training_data_set] = [1]
data[:randomize_training_data_set] = [0]

data[:fraction_test_data_set] = [1]
data[:randomize_test_data_set] = [0]

data[:basename_for_output_files] = ["HfO2FPOD"]

data[:compute_pod_descriptors] = [0]
data[:compute_pod_coefficients] = [1]


# Run fitting experiments ######################################################

run(`mkdir -p $experiments_path`)
i = 1
ks_pars, ks_dat = keys(params), keys(data)
for vs_dat in product(values(data)...)
    curr_data = OrderedDict(ks_dat .=> vs_dat)
    for vs_pars in zip(values(params)...)
        curr_params = OrderedDict(ks_pars .=> vs_pars)
    
        run(`mkdir -p $experiments_path/$i`)
        
        # data.pod
        open("$experiments_path/$i/data.pod", "w") do io
            [println(io, "$k $v") for (k, v) in curr_data]
        end

        # params.pod
        open("$experiments_path/$i/params.pod", "w") do io
            [println(io, "$k $v") for (k, v) in curr_params]
        end
        
        # fit.pod
        open("$experiments_path/$i/fit.pod", "w") do io
            println(io, "fitpod params.pod data.pod")
        end
        
        # fit pod using lammps
        if parallel
            @async run(Cmd(`nohup mpirun -n 1 $lammps_path -in fit.pod`,
                       dir="./$experiments_path/$i"))
        else
            run(Cmd(`mpirun -n 1 $lammps_path -in fit.pod`,
                dir="./$experiments_path/$i"))
        end
        
        global i += 1

    end

end



