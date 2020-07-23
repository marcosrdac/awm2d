using InteractiveUtils
include("./reff.jl")
include("./parameters.jl")

filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
io = open(filename, "w")
write(io, 2, size(v)..., v)
close(io)
