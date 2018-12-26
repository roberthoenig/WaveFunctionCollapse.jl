using FileIO
using StatsBase

input = []  # Array{RGB{Normed{UInt8,8}},2}
patternToId = Dict()
idToPattern = Dict()
patternCount = Dict()  # patternId => Int
patternAdjacency = Dict()  # Set([(direction, patternId1, patternId2), ...])

const directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]

function load(;
    filename,
    patternsize=2,
    width=16,  # Width of output image in number of patterns.
    height=16, # Height of output image in number of patterns.
    periodicInput=false,
    periodicOutput=false)
    
    global input = []
    global patternToId = Dict()
    global idToPattern = Dict()
    global patternCount = Dict()
    global patternAdjacency = Dict()

    begin  # Generate patterns.
    global input = FileIO.load(filename)
    (input_height, input_width) = size(input)
    # TODO: Account for periodicInput.
    for col in 1:input_width+1-patternsize, row in 1:input_height+1-patternsize
        pattern = input[row:row+patternsize-1, col:col+patternsize-1]
        if pattern in keys(patternToId)
            patternCount[patternToId[pattern]] += 1
        else
            newid = length(idToPattern) + 1
            idToPattern[newid] = pattern
            patternToId[pattern] = newid
            patternCount[newid] = 1
        end
    end
    end

    @assert all([patternToId[idToPattern[i]] == i for i in 1:length(patternToId)])

    begin  # Initialize patternAdjacency.
    for direction in directions, p1 in keys(idToPattern), p2 in keys(idToPattern)
        (dx, dy) = direction

        lap1 = idToPattern[p1][(dy==1 ? 2 : 1):(dy==-1 ? end-1 : end), (dx==1 ? 2 : 1):(dx==-1 ? end-1 : end)]
        lap2 = idToPattern[p2][(dy==-1 ? 2 : 1):(dy==1 ? end-1 : end), (dx==-1 ? 2 : 1):(dx==1 ? end-1 : end)]
        overlap = lap1 == lap2
        if overlap
            push!(patternAdjacency, (direction, p1, p2))
        end
    end
    end

    begin  # Generate output.
        field_patterns = [Set(keys(idToPattern)) for i in 1:height, j in 1:width]
        
        while !all(x -> size(x) == 1, field_patterns)
            if any(x -> size(x) == 0, field_patterns)
                error("Impossible to complete generation. Restart the algorithm.")
            end
            # Find field with lowest entropy in field_patterns.
            weight_sum = sum(field_patterns) do valid_patterns
                sum(id -> patternCount[id], valid_patterns)
            end
            field_entropies = map(field_patterns) do valid_patterns
                # Ignore already collapsed fields.
                if length(valid_patterns) <= 1
                    return 10e9
                end
                randn()*1e-6 - sum(valid_patterns) do id
                    weight = patternCount[id]
                    weight * log(weight/weight_sum)
                end
            end
            (_, idx_min) = findmin(field_entropies)
            field = collect(field_patterns[idx_min])
            choice = sample(field, getindex.(Ref(patternCount), field))
            field_patterns[idx_min] = Set([choice])

            # Enforce new constraint in all other fields.
            field_idxs = []
            function neighbors(idx)
                (idx+CartesianIndex(direction) for direction in directions if checkbounds(Bool, field_patterns, idx+CartesianIndex(direction)))
            end
            for neighbor in neighbors(idx)
                push!(field_idx, neighbor)
            end
            while !empty(field_idxs)
                field_idx = pop!(field_idxs)
                constrained_valid_patterns = Set([])
                for id in field_patterns[field_idx]
                    if all(neighbors(field_idx)) do constraint_idx
                        any(Base.Iterators.product(field_patterns[field_idx], field_patterns[constraint_idx])) do (field_id, constraint_id)
                            (direction, field_id, constraint_id) in patternAdjacency
                        end
                    end
                        push!(constrained_valid_patterns, id)
                    end
                end
                if length(constrained_valid_patterns) != length(field_patterns[field_idx])
                    for neighbor in neighbors(field_idx)
                        push!(field_idxs, idx)
                    end
                end
                field_patterns[field_idx] = constrained_valid_patterns
            end
        end
    end
    output_patterns = map(x -> idToPattern[first(x)], field_patterns)
    output = vcat([hcat(output_patterns[row, :]...) for row in 1:height]...)
    output
end