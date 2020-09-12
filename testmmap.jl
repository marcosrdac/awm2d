include("./propagate_module.jl")
using .Propagate
using Mmap

io = open("test.bin", "w+")

a = discarray(io, "w+",  Float64, (5,3))
a .= 1
@show a
close(io)

io = open("test.bin", "r+")
@show discarray(io, "r+")
close(io)
