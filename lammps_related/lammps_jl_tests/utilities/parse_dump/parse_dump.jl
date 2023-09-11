

"""
Question: Why doesn't the eachline iterator need the state to iterate past the first one?
"""

file_itr = eachline("./dump_single.custom")

configs = []
while !isempty(file_itr)
    config_dict = Dict()

    (line,st) = iterate(file_itr)
    @assert line == "ITEM: TIMESTEP"
    tstep = parse(Int32, iterate(file_itr)[1])

    (line,st) = iterate(file_itr)
    @assert line == "ITEM: NUMBER OF ATOMS"
    num_atoms = parse(Int32, iterate(file_itr)[1])

    (line,st) = iterate(file_itr)
    line_toks = split(line,r"\s+")
    @assert line_toks[2] == "BOX"
    boundary_style = line_toks[4:6] # type SubString{String} which maybe isn't what I want
    
    lv_bounds = []
    for _ in 1:3
        (line,st) = iterate(file_itr)
        lohi = parse.(Float64,split(line,r"\s+"))
        push!(lv_bounds,lohi)
    end

    (line,st) = iterate(file_itr)
    # This is the most fragile/specific bit. Needs pos,forces,velocities
    @assert line == "ITEM: ATOMS id type x y z fx fy fz vx vy vz"
    
    pos = []
    types = []
    forces = []
    vels = []
    last_atomid = 0 
    for _ in 1:num_atoms
        line_toks = split(strip(iterate(file_itr)[1]), r"\s+")
        @show line_toks

        # Want to ensure that they are in order (not always the case!)
        atomid = parse(Int32,line_toks[1])
        @assert atomid = last_atomid + 1
        last_atomid = atomid

    end
    
    break 



end
#    
#    
#
#
#header = true
#for line in file_itr
#  if header
#    
#
#
#for (i,line) in enumerate(file_itr)
#    @show line
#
#    if i == 5
#      for _ in 1:3
#
#        println("HEY!")
#        @show iterate(file_itr)
#      end
#    end
#
#    if i > 10
#      break
#    end 
#end
##function parse_dump(dump_fname::String)
#  
