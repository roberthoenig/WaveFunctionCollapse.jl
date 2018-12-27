using FileIO
using Random
using StatsBase

input = []  # Array{RGB{Float32},2}
patternToId = Dict()
idToPattern = Dict()
patternCount = Dict()  # patternId => Int
patternAdjacency = Set()  # Set([(direction, patternId1, patternId2), ...])

const directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]

function get_image(field_patterns; average_superpositions=false)
    output_patterns = map(field_patterns) do pattern_ids
        if average_superpositions
            error("Superposition averaging not implemented!")
            # sum(getindex.(Ref(idToPattern), collect(pattern_ids))) ./ length(pattern_ids)
        else
            pattern = length(pattern_ids) == 1 ? idToPattern[first(pattern_ids)] : zero(idToPattern[first(pattern_ids)])
            fill(pattern[1], 1, 1)
        end
    end
    vcat([hcat(output_patterns[row, :]...) for row in 1:size(output_patterns)[1]]...)
end

function print_field_patterns(field_patterns)
    println("field_patterns:")
    for idx in CartesianIndices(field_patterns)
        println("$idx: $(sort(collect(field_patterns[idx])))")
    end
end

function load(;
    filename,
    patternsize=2,
    width=16,  # Width of output image in number of patterns.
    height=16, # Height of output image in number of patterns.
    periodicInput=false,
    periodicOutput=false,
    seed=0)
    
    global input = []
    global patternToId = Dict()
    global idToPattern = Dict()
    global patternCount = Dict()
    global patternAdjacency = Set()

    if seed != 0
        Random.seed!(seed)
    else
        newseed = rand(1:1_000_000)
        Random.seed!(newseed)
        println("Seed: $newseed")
    end

    begin  # Generate patterns.
    global input = FileIO.load(filename)
    input = convert(Array{RGB{Float32}, 2}, input)
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
        (dy, dx) = direction
        lap1 = idToPattern[p1][(dy==1 ? 2 : 1):(dy==-1 ? end-1 : end), (dx==1 ? 2 : 1):(dx==-1 ? end-1 : end)]
        lap2 = idToPattern[p2][(dy==-1 ? 2 : 1):(dy==1 ? end-1 : end), (dx==-1 ? 2 : 1):(dx==1 ? end-1 : end)]
        overlap = lap1 == lap2
        if overlap
            push!(patternAdjacency, (direction, p1, p2))
        end
    end
    end
    # for direction in directions, p1 in keys(idToPattern), p2 in keys(idToPattern)
    #     (dy, dx) = direction
    #     N = patternsize
    #     xmin = dx < 0 ? 0 : dx
    #     xmax = dx < 0 ? dx + N : N
    #     ymin = dy < 0 ? 0 : dy
    #     ymax = dy < 0 ? dy + N : N
    #     overlap = true
    #     for y in ymin:ymax
    #         for x in xmin:xmax
    #             if p1[y, x] != p2[y-dy, x-dx]
    #                 overlap = false
    #             end
    #         end
    #     end
    #     if overlap
    #         push!(patternAdjacency, (direction, p1, p2))
    #     end
    # end


    begin  # Generate output.
        field_patterns = [Set(keys(idToPattern)) for i in 1:height, j in 1:width]
        while !all(x -> length(x) == 1, field_patterns)
            # print_field_patterns(field_patterns)
            # flush(stdout)
            if any(x -> length(x) == 0, field_patterns)
                error("Impossible to complete generation. Restart the algorithm.")
            end
            # display(get_image(field_patterns, average_superpositions=false))
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
            weights = Weights(getindex.(Ref(patternCount), field))
            # @show weights
            choice = sample(field, weights)
            field_patterns[idx_min] = Set([choice])
            # println("Choosing pattern $choice")
            # println("finished field selection")
            flush(stdout)
            # Enforce new constraint in all other fields.
            field_idxs = []
            function neighbors(idx)
                (idx+CartesianIndex(direction) for direction in directions if checkbounds(Bool, field_patterns, idx+CartesianIndex(direction)))
            end
            for neighbor in neighbors(idx_min)
                push!(field_idxs, neighbor)
            end
            while !isempty(field_idxs)
                field_idx = pop!(field_idxs)
                constrained_valid_patterns = Set([])
                for pattern_id in field_patterns[field_idx]
                    if all(directions) do direction
                        constraint_idx = field_idx+CartesianIndex(direction)
                        if !checkbounds(Bool, field_patterns, constraint_idx) return true end
                        # println("comparing field $field_idx with $constraint_idx")
                        any(field_patterns[constraint_idx]) do constraint_pattern_id
                            (direction, pattern_id, constraint_pattern_id) in patternAdjacency
                        end
                    end
                        push!(constrained_valid_patterns, pattern_id)
                    end
                end
                if length(constrained_valid_patterns) != length(field_patterns[field_idx])
                    for neighbor in neighbors(field_idx)
                        push!(field_idxs, neighbor)
                    end
                end
                field_patterns[field_idx] = constrained_valid_patterns
                # println("constraint enforcing")
                # print_field_patterns(field_patterns)
            end
        end
    end
    # print_field_patterns(field_patterns)
    println("finished")
    get_image(field_patterns)
end