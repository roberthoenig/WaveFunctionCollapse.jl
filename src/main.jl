include("WaveFunctionCollapse.jl")

using ArgParse
using FileIO
using .WaveFunctionCollapse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "filename"
            required = true
        "--patternsize"
            arg_type = Int
            default = 3
        "--width"
            arg_type = Int
            default = 32
        "--height"
            arg_type = Int
            default = 32
        "--periodic_input"
            action = :store_true
        "--periodic_output"
            action = :store_true
        "--mirror_input_horizontally"
            action = :store_true
        "--mirror_input_vertically"
            action = :store_true
        "--rotate_input_clockwise"
            action = :store_true
        "--rotate_input_anticlockwise"
            action = :store_true
        "--seed"
            arg_type = Int
            default = 0
        "--save_to_gif"
            action = :store_true
        "--fast"
            action = :store_true
    end

    return parse_args(s, as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    filename = parsed_args[:filename]
    delete!(parsed_args, :filename)
    img = generate(filename; parsed_args...)
    FileIO.save("output.png", img)
    println("Done! Saved the output image as output.png.")
end

main()