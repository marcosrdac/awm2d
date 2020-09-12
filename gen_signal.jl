include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 1 # Hz

signature = rickerwave(ν, Δt)

todiscarray(sourcesignaturefile, signature)
