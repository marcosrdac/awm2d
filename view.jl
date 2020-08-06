include("./propagate_module.jl")
include("./parameters.jl")

using PyPlot
using .Propagate

#filename = "data/seis.bin"
filename, nt = P_file, 300


if occursin("P", filename)
    P = discarray(filename, "r")
    plt.imshow(P[:,:,nt])
elseif occursin("seis", filename)
    seis = discarray(filename, "r")
    plt.imshow(seis; aspect="auto", vmin=-.006, vmax=.006)
end

println(seis)
#plt.show()
