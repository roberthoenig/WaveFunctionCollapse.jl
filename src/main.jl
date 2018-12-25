include("wavefunctioncollapse.jl")

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-n", "--patternsize"
            help = "size (width and height) of the patterns"
            arg_type = Int
            default = 2
        "--width"
            help = "width of the output image"
            arg_type = Int
            default = 64
        "--height"
            help = "height of the output image"
            arg_type = Int
            default = 64
        "--periodicInput"
            help = "whether the input image should be treated as periodic"
            action = :store_true
        "--periodicOutput"
            help = "whether the output image should be generated to be periodic"
            action = :store_true
        "filename"
            help = "path to the input image"
            required = true
    end

    return parse_args(s, as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    load(; parsed_args...)
end

main()