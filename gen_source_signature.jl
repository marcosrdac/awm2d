include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

ν = 10 # Hz

signature = rickerwave(ν, Δt)

todiscarray(sourcesignaturefile, signature)
