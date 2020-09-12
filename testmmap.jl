include("./propagate_module.jl")
using .Propagate

io = open("test.bin", "w+")
# io = "test.bin"

# @show discarray(io, "w+",  Float64, (5,3)) .= rand()
@show todiscarray(io, rand(Float64, (3,3)))
println()
@show discarray(io, "r+"; pos=0)
