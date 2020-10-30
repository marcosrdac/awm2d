include("src/acoustics2d.jl")
include("src/discarrays.jl")
include("parameters.jl")
using .Acoustics2D
using .Discarrays


ν = .5 # Hz


signature = rickerwave(ν, Δt)

todiscarray(sourcesignaturefile, signature)
