include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 10 # Hz
source_signature = todiscarray(source_signature_file, rickerwave(ν, Δt))
