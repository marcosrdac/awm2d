include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 60 # Hz
source_signature = todiscarray(source_signature_file, rickerwave(ν, Δt))
