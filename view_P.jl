using PyPlot
using Mmap: mmap


function open_mmap(filename::String, mode::String="r")
    io = open(filename, mode)
    _ndims = read(io, Int64)
    _dims = Tuple(read(io, Int64) for i in 1:_ndims)
    P = mmap(io, Array{Float64, _ndims}, _dims)
    close(io)
    return(P)
end


#P = open_mmap("P.bin", "r")
#plt.imshow(P[:,:,1500])


seis = open_mmap("seis.bin", "r")
plt.imshow(seis; aspect="auto", vmin=-.006, vmax=.006)

plt.show()
