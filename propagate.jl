using InteractiveUtils
include("./reff.jl")


# mesh and solving parameters
begin
    h  = 1.0 # km
    Δt = .001 # s
    NX = 321
    NZ = 321
    NT = 2700
    grid = FDM_Grid(h, Δt, NZ, NX, NT)
end


# three layered model parameters
begin
    H1 = (2*NZ)÷8
    H2 = (3*NZ)÷8
    V1 = 3. # km/s
    V2 = 4. # km/s
    V3 = 5. # km/s

    # reflectors position
    z1 = H1
    z2 = z1+H2

    # actually defining velocity field
    v = Array{Float64}(undef, (NZ,NX))
    v[   1:z1,  1:end] .= V1
    v[z1+1:z2,  1:end] .= V2
    v[z2+1:end, 1:end] .= V3
end


# signal parameters
begin
    ν = 3 # Hz
    array = "split"
    signature = rickerwave(ν, Δt)

    if array === "split"      position = CartesianIndex(1, NX÷2)
    elseif array === "endon"  position = CartesianIndex(1, 1)
    elseif array === "center" position = CartesianIndex(NZ÷2, NX÷2)
    end
    signal = Signal(signature, position)
end


# starting pressure field
begin
    P0 = zero(v)
    #P0[1,1] = 1      # endon test
    #P0[1, NX÷2] = 1  # split test
end


## padding arrays, allocating 2 more times for pressure field
#begin
#    _v = pad_extremes(v, TAPER)
#    _P = pad_zeros_add_axes(P0, TAPER+1, 3)
#end


#@time propagate_absorb(grid, _P, _v,  signal)
#@time P = propagate_save(grid, P0, v,  signal; filename="P.bin")
#@time propagate_save(grid, P0, v, signal; filename="P.bin")
@time save_seis(grid, P0, v, signal; filename="seis.bin")
#@time propagate_absorb(grid, _P, _v,  signal)


#using BenchmarkTools
### runs 6-7 times, be careful
#print("Propagate absorb ")
#@btime propagate_absorb($grid, $_P, $_v,  $signal)  # 10.986
#print("Propagate pure ")
#@btime propagate_pure($grid, $_P, $_v,  $signal)




#using PyPlot
#
##plt.imshow(P[1+TAPER:end-TAPER, 1+TAPER:end-TAPER, 1])
##plt.imshow(P[:, :, 1]; vmin=-.2, vmax=.2)
##plt.imshow(P[:, :, 1])
##plt.imshow(P[:, :])
#
#plt.colorbar()
#plt.show()
