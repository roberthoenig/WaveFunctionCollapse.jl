using FileIO
using Random
using StatsBase
using FFTViews
using Images

export generate

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
PeriodicArray = FFTView

function getImage(field_patterns, id_to_pattern; average_superpositions=false)
    output_patterns = map(field_patterns) do pattern_ids
        if average_superpositions
            pattern = sum(getindex.(Ref(id_to_pattern), collect(pattern_ids)))[end, 1] / length(pattern_ids)
        else
            pattern = length(pattern_ids) == 1 ? id_to_pattern[first(pattern_ids)] : zero(id_to_pattern[first(pattern_ids)])
            fill(pattern[end, 1], 1, 1)
        end
    end
    vcat([hcat(output_patterns[row, :]...) for row in 1:size(output_patterns)[1]]...)
end

const directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
const opposite = Dict((1,0) => (-1,0), (0,1) => (0, -1), (-1, 0) => (1, 0), (0, -1) => (0, 1))

function getNeighbors(idx, field_patterns)
    ((idx+CartesianIndex(direction), direction) for direction in directions
        if checkbounds(Bool, field_patterns, idx+CartesianIndex(direction)))
end
function disallowPattern(field_idx, pattern_id, field_patterns, pattern_adjacency, patterns_allowed)
    for (neighbor, direction) in getNeighbors(field_idx, field_patterns),
            adj_pattern_id in pattern_adjacency[direction][pattern_id]
        patterns_allowed[opposite[direction]][Tuple(neighbor)..., adj_pattern_id] -= 1
    end
end

function expandFrontier(field_idx, field_patterns, frontier)
    offset = 3
    min_row = field_idx[1]-offset
    max_row = field_idx[1]+offset
    min_col = field_idx[2]-offset
    max_col = field_idx[2]+offset
    surrounding_fields = [CartesianIndex(row, col) for (row, col) in Iterators.product(min_row:max_row, min_col:max_col)]
    surrounding_fields = filter(f -> checkbounds(Bool, field_patterns, f), surrounding_fields)
    for field in surrounding_fields
        if length(field_patterns[field]) > 1
            push!(frontier, field)
        end
    end
    delete!(frontier, field_idx)
end

function collapseField(collapse_field_idx, chosen_pattern_id, field_patterns, pattern_adjacency, patterns_allowed, fast, frontier)
    #  Constrain patterns_allowed with the new information.
    if fast
        expandFrontier(collapse_field_idx, field_patterns, frontier)
    end
    field = collect(field_patterns[collapse_field_idx])
    for pattern_id in field
        if pattern_id == chosen_pattern_id continue end
        disallowPattern(collapse_field_idx, pattern_id, field_patterns, pattern_adjacency, patterns_allowed)
    end
    field_patterns[collapse_field_idx] = Set([chosen_pattern_id])
    # Enforce new constraint in all other fields.
    field_idxs = [neighbor for (neighbor, _) in getNeighbors(collapse_field_idx, field_patterns)]
    while !isempty(field_idxs)
        field_idx = pop!(field_idxs)
        constrained_valid_patterns = Set([])
        for pattern_id in field_patterns[field_idx]
            if all(getNeighbors(field_idx, field_patterns)) do (_, direction)
                    patterns_allowed[direction][field_idx[1], field_idx[2], pattern_id] > 0
                end
                push!(constrained_valid_patterns, pattern_id)
            else
                disallowPattern(field_idx, pattern_id, field_patterns, pattern_adjacency, patterns_allowed)
            end
        end
        if length(constrained_valid_patterns) != length(field_patterns[field_idx])
            for (neighbor, _) in getNeighbors(field_idx, field_patterns)
                if (!fast) || (neighbor in frontier)
                    push!(field_idxs, neighbor)
                end
            end
        end
        field_patterns[field_idx] = constrained_valid_patterns
        if fast && length(field_patterns[field_idx]) == 1
            expandFrontier(field_idx, field_patterns, frontier)
        end
    end
end

"""
    generate(filename::String
        [; patternsize=2,
           width=16,
           height=16,
           periodic_input=false,
           periodic_output=false,
           mirror_input_horizontally=false,
           mirror_input_vertically=false,
           rotate_input_clockwise=false,
           rotate_input_anticlockwise=false,
           ground=nothing,
           seed=0,
           save_to_gif=false,
           fast=false]
    )

Apply the WaveFunctionCollapse algorithm to generate an image of dimension `width`*`height`
based on the input image at location `filename`. The algorithm splits the input image into
many tiles and then constructs the output image by places similar tiles next to each other.
`patternsize` is the width and height of a pattern. If `periodic_input` is true, the algorithm
will consider patterns in the input image that wrap around vertically and horizontally. If
`periodic_output` is true, the output image can be seamlessly repeated vertically and horizontally.
If the `mirror_input_` and `rotate_input_` parameters are true, the algorithm generates patterns
also from the mirrored and rotated input image. If `ground` is a tuple (y, x), the left bottom field
of the output is set to the input pattern with the upper left corner at (y, x). In practice, this
can enforce regularity on images with ground rows, like `samples/Flowers.png`. If `seed` != 0,
calling the algorithm with the same seed will generate the same output image, provided that the other
parameters are identical. If `save_to_gif` is true, the generation process gets animated and saved
as "output.gif".

If `fast` is true, an experimental version of the algorithm is used. In practice, it is 2-4 times as fast,
but may lead to more failed generation attempts when generating big outputs. The experimental algorithm is
almost like the normal algorithm, but it doesn't propagate pattern constraints through all fields. Instead,
it limits the propagation to a frontier around already collapsed fields. That frontier is currently 3 fields
thick, but that value is subject to experimentation.

    generate(input::Array{<:ColorTypes.RGB,2}
    [; kwargs...]
    )

Like `generate(filename [; kwargs...])`, but you can directly pass an array for the input image.

# Examples
```
generate(
    filename="samples/Flowers.png",
    width=20,
    height=20,
    patternsize=3,
    periodic_output=true,
    mirror_input_horizontally=true,
    seed=726259)
```
"""
function generate(input::Array{<:ColorTypes.RGB,2};
    patternsize=2,
    width=16,
    height=16,
    periodic_input=false,
    periodic_output=false,
    mirror_input_horizontally=false,
    mirror_input_vertically=false,
    rotate_input_clockwise=false,
    rotate_input_anticlockwise=false,
    ground=nothing,
    seed=0,
    save_to_gif=false,
    fast=false,
)
    generate(input, patternsize, width, height, periodic_input,
        periodic_output, mirror_input_horizontally, mirror_input_vertically,
        rotate_input_clockwise, rotate_input_anticlockwise, ground, seed,
        save_to_gif, fast
    )
end

generate(filename::String; kwargs...) = generate(FileIO.load(filename); kwargs...)

function generate(
    input::Array{<:ColorTypes.RGB,2}, patternsize, width, height, periodic_input, periodic_output,
    mirror_input_horizontally, mirror_input_vertically, rotate_input_clockwise,
    rotate_input_anticlockwise, ground, seed, save_to_gif, fast,
)
    input = convert(Array{RGB{Float32}, 2}, input)
    pattern_to_id = Dict()
    id_to_pattern = Dict()
    pattern_count = Dict{Int64, Int64}()
    images = []
    frontier = Set{CartesianIndex{2}}()

     # direction => pattern_id => Set([pattern_ids...])
    pattern_adjacency = Dict{Tuple{Int64,Int64},Dict{Int64,Set{Int64}}}()

    # direction => [width x height x patterns]. Pattern i at (x, y) overlaps with patterns
    # from the field at (x, y)+direction patterns_allowed[direction][x, y, i] times.
    patterns_allowed = Dict{Tuple{Int64,Int64},Array{Int64,3}}() 

    if seed != 0
        Random.seed!(seed)
    else
        newseed = rand(1:1_000_000)
        Random.seed!(newseed)
        println("Seed: $newseed")
    end
    bound = patternsize-1
    if periodic_input
        input = PeriodicArray(input)
        bound = 0
    end
    ground_pattern_id = nothing
    for col in 1:size(input)[2]-bound, row in 1:size(input)[1]-bound
        pattern_variations = [input[row:row+patternsize-1, col:col+patternsize-1]]
        if mirror_input_horizontally
            append!(pattern_variations, map(p -> p[1:end, end:-1:1], pattern_variations))
        end
        if mirror_input_vertically
            append!(pattern_variations, map(p -> p[end:-1:1, 1:end], pattern_variations))
        end
        if rotate_input_clockwise
            append!(pattern_variations, map(rotr90, pattern_variations))
        end
        if rotate_input_anticlockwise
            append!(pattern_variations, map(rotl90, pattern_variations))
        end
        for pattern in pattern_variations
            if pattern in keys(pattern_to_id)
                pattern_count[pattern_to_id[pattern]] += 1
            else
                newid = length(id_to_pattern) + 1
                id_to_pattern[newid] = pattern
                pattern_to_id[pattern] = newid
                pattern_count[newid] = 1
            end
        end
        if (row, col) == ground
            ground_pattern_id = pattern_to_id[pattern_variations[1]]
            println("Ground:"); flush(stdout); display(id_to_pattern[ground_pattern_id])
        end
    end

    for direction in directions, p1 in keys(id_to_pattern), p2 in keys(id_to_pattern)
        (dy, dx) = direction
        lap1 = id_to_pattern[p1][(dy==1 ? 2 : 1):(dy==-1 ? end-1 : end), (dx==1 ? 2 : 1):(dx==-1 ? end-1 : end)]
        lap2 = id_to_pattern[p2][(dy==-1 ? 2 : 1):(dy==1 ? end-1 : end), (dx==-1 ? 2 : 1):(dx==1 ? end-1 : end)]
        overlap = lap1 == lap2
        if overlap
            # TODO: Properly initialize pattern_adjacency
            push!(get!(get!(pattern_adjacency, direction, Dict()), p1, Set()), p2)
        end
    end
    for direction in directions
        # WARNING: Requires that pattern IDs go from 1 to n.
        patterns = map(1:length(id_to_pattern)) do id
            length(get!(pattern_adjacency[direction], id, Set()))
        end
        patterns_allowed[direction] = [p for _ in 1:height, _ in 1:width, p in patterns]
        if periodic_output
            patterns_allowed[direction] = PeriodicArray(patterns_allowed[direction])
        end
    end

    field_patterns = [Set(keys(id_to_pattern)) for i in 1:height, j in 1:width]
    if periodic_output
        field_patterns = PeriodicArray(field_patterns)
    end
    if ground_pattern_id != nothing
        collapseField(CartesianIndex(height, 1), ground_pattern_id, field_patterns, pattern_adjacency, patterns_allowed, fast, frontier)
        if any(x -> length(x) == 0, field_patterns)
            error("Impossible to initialize bottom row with ground. Try another ground.")
        end
    end
    while !all(x -> length(x) == 1, field_patterns)
        # Find field with lowest entropy in field_patterns.
        weight_sum = sum(field_patterns) do valid_patterns
            sum(id -> pattern_count[id], valid_patterns)
        end
        field_entropies = map(field_patterns) do valid_patterns
            # Ignore already collapsed fields.
            if length(valid_patterns) <= 1
                return 10e9
            end
            randn()*1e-6 - sum(valid_patterns) do id
                weight = pattern_count[id]
                weight * log(weight/weight_sum)
            end
        end
        (_, idx_min) = findmin(field_entropies)
        field = collect(field_patterns[idx_min])
        choice = sample(field, Weights(getindex.(Ref(pattern_count), field)))
        collapseField(idx_min, choice, field_patterns, pattern_adjacency, patterns_allowed, fast, frontier)
        if any(x -> length(x) == 0, field_patterns)
            error("Impossible to complete generation. Restart the algorithm.")
        end
        if save_to_gif
            push!(images, getImage(field_patterns, id_to_pattern, average_superpositions=true))
        end
    end
    if save_to_gif
        save("output.gif", cat(images..., dims=[3]), fps=10)
    end
    getImage(field_patterns, id_to_pattern)
end