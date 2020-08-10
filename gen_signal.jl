include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 6 # Hz
source_signature = todiscarray(source_signature_file, rickerwave(ν, Δt))
