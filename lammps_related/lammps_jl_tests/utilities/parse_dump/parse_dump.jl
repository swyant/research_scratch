#using Printf

"""
Question: Why doesn't the eachline iterator need the state to iterate past the first one?
"""

function fragile_parse_dump(fname::String)
    file_itr = eachline(fname)
    
    configs = []
    while !isempty(file_itr)
        config_dict = Dict()
    
        (line,st) = iterate(file_itr)
        @assert line == "ITEM: TIMESTEP"
        tstep = parse(Int32, iterate(file_itr)[1])
        config_dict["timestep"] = tstep
    
        (line,st) = iterate(file_itr)
        @assert line == "ITEM: NUMBER OF ATOMS"
        num_atoms = parse(Int32, iterate(file_itr)[1])
        config_dict["num_atoms"] = num_atoms
    
        (line,st) = iterate(file_itr)
        line_toks = split(line,r"\s+")
        @assert line_toks[2] == "BOX"
        boundary_style = line_toks[4:6] # type SubString{String} which maybe isn't what I want
        config_dict["boundary"] = copy(boundary_style)
        
        lv_bounds = []
        for _ in 1:3
            (line,st) = iterate(file_itr)
            lohi = parse.(Float64,split(line,r"\s+"))
            push!(lv_bounds,lohi)
        end
    
        config_dict["lv_bounds"] = deepcopy(lv_bounds)
    
        (line,st) = iterate(file_itr)
        # This is the most fragile/specific bit. Needs pos,forces,velocities
        @assert line == "ITEM: ATOMS id type x y z fx fy fz vx vy vz"
        
        cfg_poss = []
        cfg_types = Int32[]
        cfg_forces = []
        cfg_vels = []
        last_atomid = 0 
        for _ in 1:num_atoms
            line_toks = split(strip(iterate(file_itr)[1]), r"\s+")
            #@show line_toks
    
            # Want to ensure that they are in order (not always the case!)
            atomid = parse(Int32,line_toks[1])
            @assert atomid == last_atomid + 1
            last_atomid = atomid
            
            cfg_type = parse(Int32, line_toks[2])
            push!(cfg_types,cfg_type)
    
            cfg_pos = parse.(Float64, line_toks[3:5])
            #@printf "%32.27f\n" cfg_pos[1]
            push!(cfg_poss,cfg_pos)
    
            cfg_force = parse.(Float64,line_toks[6:8])
            #@printf "%32.27f\n" cfg_force[1]
            push!(cfg_forces,cfg_force)
    
            cfg_vel  = parse.(Float64,line_toks[9:11])
            #@printf "%32.27f\n" cfg_vel[1]
            push!(cfg_vels, cfg_vel)
        end
    
        config_dict["positions"]  = reduce(hcat,cfg_poss)
        config_dict["types"]      = cfg_types
        config_dict["forces"]     = reduce(hcat,cfg_forces)
        config_dict["velocities"] = reduce(hcat, cfg_vels)
    
        
        push!(configs,config_dict)
   
    end   
    configs
end 

#test_configs = fragile_parse_dump("./dump_single.custom")
