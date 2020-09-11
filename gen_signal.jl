include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 10 # Hz
source_signature = todiscarray(sourcesignaturefile, rickerwave(ν, Δt))
