using InteractiveUtils
include("./reff.jl")
include("./parameters.jl")


# mesh and solving parameters
begin
    # h  = 1.0 # km
    # Δt = .001 # s
    # NX = 321
    # NZ = 321
    # NT = 3000
    grid = FDM_Grid(h, Δt, NZ, NX, NT)
end


# three layered model parameters
# begin
    # H1 = (1*NZ)÷8
    # H2 = (1*NZ)÷8
    # V1 = 3. # km/s
    # V2 = 4. # km/s
    # V3 = 5. # km/s

    # # reflectors position
    # z1 = H1
    # z2 = z1+H2

    # # actually defining velocity field
    # v = Array{Float64}(undef, (NZ,NX))
    # v[   1:z1,  1:end] .= V1
    # v[z1+1:z2,  1:end] .= V2
    # v[z2+1:end, 1:end] .= V3
# end


# signal parameters
begin
    ν = 6 # Hz
    array = "split"
    signature = rickerwave(ν, Δt)

    if array === "split"      position = CartesianIndex(1, NX÷2+1)
    elseif array === "endon"  position = CartesianIndex(1, 1)
    elseif array === "center" position = CartesianIndex(NZ÷2+1, NX÷2+1)
    end
    signal = Signal1D(signature, position)
end


# starting pressure field
begin
    P0 = zero(v)
end


@time P = propagate(grid, P0, v, signal;
                   save=true,
                   filename="/mnt/hdd/home/tmp/awp_data/P.bin",
                   only_seis=false)

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
#
###plt.imshow(P[1+TAPER:end-TAPER, 1+TAPER:end-TAPER, 1])
##plt.imshow(P[:, :, 1]; vmin=-.2, vmax=.2)
#plt.imshow(P[:, :, 1])
##plt.imshow(P[:, :])
#
#plt.colorbar()
#plt.show()
