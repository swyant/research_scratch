using DelimitedFiles
using JLD2 

dir = "./2body_rerun_monoclinic_hfo2"
out = "monoclinic_hfo2_lammps_state.jld2"

function read_ragged_array(filename, t)
    result = Vector{Vector{t}}()
    
    open(filename) do file
        for line in eachline(file)
            # Split the line by whitespace and convert each value to Float64
            # Skip empty strings that might result from extra spaces
            row = [parse(t, val) for val in split(line) if !isempty(val)]
            push!(result, row)
        end
    end    
    return result
end


x_file = "$(dir)/x.txt"
ilist_file = "$(dir)/ilist.txt"
numneigh_file = "$(dir)/numneigh.txt"
firstneigh_file = "$(dir)/firstneigh.txt"
type_file = "$(dir)/type.txt"

x =readdlm(x_file)
ilist = vec(readdlm(ilist_file, Int64))
numneigh = vec(readdlm(numneigh_file, Int64))
type = vec(readdlm(type_file, Int64))
firstneigh = read_ragged_array(firstneigh_file,Int64)
inum = length(type) 
@assert inum == length(ilist) == length(numneigh)

save(out, Dict("x" => x,
               "type" => type, 
               "inum" => inum,
               "ilist" => ilist,
               "numneigh" => numneigh,
               "firstneigh" => firstneigh))
