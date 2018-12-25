# using Images
using FileIO

input = []
patternToId = Dict()
idToPattern = Dict()
patternCount = Dict()
patternAdjacency = Dict()

const directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]


function load(;
    filename,
    patternsize=2,
    width=64,  # Output.
    height=64, # Output.
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
            adjacency = get!(patternAdjacency, direction, Dict())
            adjacent = get!(adjacency, p1, Set())
            push!(adjacent, p2)
        end
    end
    end
end
