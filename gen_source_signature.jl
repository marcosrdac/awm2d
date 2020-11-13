include("src/acoustics2d.jl")
include("src/discarrays.jl")
include("parameters.jl")
using .Acoustics2D
using .Discarrays


ν = 15 # Hz

signature = rickerwave(ν, Δt)

todiscarray(sourcesignaturefile, signature)
