using CSV, DataFrames 

pe_data = CSV.read("./pe.dat", DataFrame, header=["tstep", "pe"], skipto=2, delim=" ")

pe_data.tstep
