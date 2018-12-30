using FileIO
using Random
using StatsBase
using FFTViews
using Images

export generate

PeriodicArray = FFTView
# Overwrite FFTViews method to use 1-based indexing.
FFTViews.indrange(i) = FFTViews.URange(first(i), last(i))
FFTViews._reindex(::Type{FFTView}, ind, i) = FFTViews.modrange(i, ind)

function Base.similar(A::AbstractArray, T::Type, shape::Tuple{FFTViews.FFTVRange,Vararg{FFTViews.FFTVRange}})
    all(x->first(x)==1, shape) || throw(BoundsError("cannot allocate FFTView with the first element of the range non-zero"))
    FFTViews.FFTView(similar(A, T, map(length, shape)))
end

function Base.similar(f::Union{Function,Type}, shape::Tuple{FFTViews.FFTVRange,Vararg{FFTViews.FFTVRange}})
    all(x->first(x)==1, shape) || throw(BoundsError("cannot allocate FFTView with the first element of the range non-zero"))
    FFTViews.FFTView(similar(f, map(length, shape)))
end

function generate(;
    filename,
    patternsize=2,
    width=16,  # Width of output image in number of patterns.
    height=16, # Height of output image in number of patterns.
    periodicInput=false,
    periodicOutput=false,
    mirrorInputHorizontally=false,
    mirrorInputVertically=false,
    rotateInputClockwise=false,
    rotateInputAnticlockwise=false,
    seed=0,
    save_to_gif=false)
    
    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    opposite = Dict((1,0) => (-1,0), (0,1) => (0, -1), (-1, 0) => (1, 0), (0, -1) => (0, 1))
    input = []  # Array{RGB{Float32},2}
    patternToId = Dict()
    idToPattern = Dict()
    patternCount = Dict()  # patternId => Int
    patternAdjacency = Dict()  # direction => pattern1 => Set([patterns...])
    patternsAllowed = Dict()  # direction => [width x height x patterns]
    #  pattern i at (x, y) overlaps with patterns from the field at (x, y)+direction patternsAllowed[direction][x, y, i] times.
    images = []

    function get_image(field_patterns; average_superpositions=false)
        output_patterns = map(field_patterns) do pattern_ids
            if average_superpositions
                pattern = sum(getindex.(Ref(idToPattern), collect(pattern_ids)))[1] / length(pattern_ids)
            else
                pattern = length(pattern_ids) == 1 ? idToPattern[first(pattern_ids)] : zero(idToPattern[first(pattern_ids)])
                fill(pattern[1], 1, 1)
            end
        end
        vcat([hcat(output_patterns[row, :]...) for row in 1:size(output_patterns)[1]]...)
    end

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
        bound = patternsize-1
        if periodicInput
            input = PeriodicArray(input)
            bound = 0
        end
        for col in 1:input_width-bound, row in 1:input_height-bound
            patternVariations = [input[row:row+patternsize-1, col:col+patternsize-1]]
            if mirrorInputHorizontally
                append!(patternVariations, map(p -> p[1:end, end:-1:1], patternVariations))
            end
            if mirrorInputVertically
                append!(patternVariations, map(p -> p[end:-1:1, 1:end], patternVariations))
            end
            if rotateInputClockwise
                append!(patternVariations, map(rotr90, patternVariations))
            end
            if rotateInputAnticlockwise
                append!(patternVariations, map(rotl90, patternVariations))
            end
            for pattern in patternVariations
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
    end

    begin  # Initialize patternAdjacency.
        for direction in directions, p1 in keys(idToPattern), p2 in keys(idToPattern)
            (dy, dx) = direction
            lap1 = idToPattern[p1][(dy==1 ? 2 : 1):(dy==-1 ? end-1 : end), (dx==1 ? 2 : 1):(dx==-1 ? end-1 : end)]
            lap2 = idToPattern[p2][(dy==-1 ? 2 : 1):(dy==1 ? end-1 : end), (dx==-1 ? 2 : 1):(dx==1 ? end-1 : end)]
            overlap = lap1 == lap2
            if overlap
                # TODO: Properly initialize patternAdjacency
                push!(get!(get!(patternAdjacency, direction, Dict()), p1, Set()), p2)
            end
        end
        for direction in directions
            # WARNING: Works only when pattern IDs go from 1 to n.
            patterns = map(1:length(idToPattern)) do id
                length(get!(patternAdjacency[direction], id, Set()))
            end
            patternsAllowed[direction] = [p for _ in 1:height, _ in 1:width, p in patterns]
            if periodicOutput
                patternsAllowed[direction] = PeriodicArray(patternsAllowed[direction])
            end
        end
    end

    begin  # Generate output.
        field_patterns = [Set(keys(idToPattern)) for i in 1:height, j in 1:width]
        if periodicOutput
            field_patterns = PeriodicArray(field_patterns)
        end
        function neighbors(idx)
            ((idx+CartesianIndex(direction), direction) for direction in directions if checkbounds(Bool, field_patterns, idx+CartesianIndex(direction)))
        end
        function disallowPattern(field_idx, pattern_id)
            for (neighbor, direction) in neighbors(field_idx)
                # TODO: This line should be unnecessary, because pattern_id should already be removed from field_idx.
                (y, x) = Tuple(field_idx)
                (neighbor_y, neighbor_x) = Tuple(neighbor)
                patternsAllowed[direction][y, x, pattern_id] = 0
                for adj_pattern_id in patternAdjacency[direction][pattern_id]
                    patternsAllowed[opposite[direction]][neighbor_y, neighbor_x, adj_pattern_id] -= 1
                end
            end
        end
        while !all(x -> length(x) == 1, field_patterns)
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
            choice = sample(field, weights)
            #  Constrain patternsAllowed with the new information.
            for pattern_id in field
                if pattern_id == choice continue end
                disallowPattern(idx_min, pattern_id)
            end
            field_patterns[idx_min] = Set([choice])
            flush(stdout)
            # Enforce new constraint in all other fields.
            field_idxs = [neighbor for (neighbor, _) in neighbors(idx_min)]
            while !isempty(field_idxs)
                field_idx = pop!(field_idxs)
                constrained_valid_patterns = Set([])
                for pattern_id in field_patterns[field_idx]
                    if all(neighbors(field_idx)) do (_, direction)
                            patternsAllowed[direction][field_idx[1], field_idx[2], pattern_id] > 0
                        end
                        push!(constrained_valid_patterns, pattern_id)
                    else
                        disallowPattern(field_idx, pattern_id)
                    end
                end
                if length(constrained_valid_patterns) != length(field_patterns[field_idx])
                    for (neighbor, _) in neighbors(field_idx)
                        push!(field_idxs, neighbor)
                    end
                end
                field_patterns[field_idx] = constrained_valid_patterns
            end
            if any(x -> length(x) == 0, field_patterns)
                error("Impossible to complete generation. Restart the algorithm.")
            end
            push!(images, get_image(field_patterns, average_superpositions=true))
        end
    end
    if save_to_gif
        save("output.gif", cat(images..., dims=[3]), fps=10)
    end
    get_image(field_patterns)
end