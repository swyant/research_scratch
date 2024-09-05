using CSV 
using DataFrames 

#gpt 3.5, small modifications
function read_segmented_files(fname::AbstractString; kwargs...)
    segs = Dict{String, DataFrame}()

    curr_seg_name = ""
    curr_seg_data = Vector{String}()

    # Read the file line by line
    for line in eachline(fname)
        if startswith(line, "# ")
            # Found a comment line indicating a new segment
            seg_name = strip(split(line, " ")[2])
            if !isempty(curr_seg_data)
                # Process the previous segment
                curr_seg_df = CSV.read(IOBuffer(join(curr_seg_data, "\n")), DataFrame; kwargs...)
                segs[curr_seg_name] = curr_seg_df
            end
            # Prepare for the new segment
            curr_seg_name = seg_name
            curr_seg_data = []
        else
            # Accumulate data for the current segment
            push!(curr_seg_data, line)
        end
    end

    # Process the last segment
    if !isempty(curr_seg_data)
        curr_seg_df = CSV.read(IOBuffer(join(curr_seg_data, "\n")), DataFrame; kwargs...)
        segs[curr_seg_name] = curr_seg_df
    end

    return segs
end

