using InteractiveUtils
include("./reff.jl")
include("./parameters.jl")

@time P = subtract_direct(P_file, direct_file, P_sub_direct_file)

#using BenchmarkTools
#@btime P = propagate_pure($grid, $P0, $v, $signal)
# save=false
# filename="tmpssd/P.bin"
# only_seis=true
# @btime P = propagate($grid, $P0, $v, $signal;
#                      save=$save,
#                      filename=$filename,
#                      only_seis=$only_seis)


#using PyPlot
###plt.imshow(P[1+TAPER:end-TAPER, 1+TAPER:end-TAPER, 1])
##plt.imshow(P[:, :, 1]; vmin=-.2, vmax=.2)
#plt.imshow(P[:, :, 1])
##plt.imshow(P[:, :])
#plt.colorbar()
#plt.show()
