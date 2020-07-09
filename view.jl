using PyPlot
using Mmap: mmap


#filename = "data/seis.bin"
filename = "data/seis.bin"
nt = 300


function open_mmap(filename::String, mode::String="r")
    io = open(filename, mode)
    _ndims = read(io, Int64)
    _dims = Tuple(read(io, Int64) for i in 1:_ndims)
    P = mmap(io, Array{Float64, _ndims}, _dims)
    close(io)
    return(P)
end


if occursin("P", filename)
    P = open_mmap(filename, "r")
    plt.imshow(P[:,:,nt])
elseif occursin("seis", filename)
    seis = open_mmap(filename, "r")
    plt.imshow(seis; aspect="auto", vmin=-.006, vmax=.006)
end


plt.show()
