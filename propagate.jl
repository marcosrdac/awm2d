using InteractiveUtils
include("./reff.jl")

println("Number of threads = $(nthreads())")


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
    v0 = Array{Float64}(undef, (NZ,NX))
    v0[   1:z1,  1:end] .= V1
    v0[z1+1:z2,  1:end] .= V2
    v0[z2+1:end, 1:end] .= V3
end


# signal parameters
begin
    ν = 6 # Hz
    array = "split"
    signature = rickerwave(ν, Δt)

    if array === "split"      position = CartesianIndex(1, NX÷2)
    elseif array === "endon"  position = CartesianIndex(1, 1)
    elseif array === "center" position = CartesianIndex(NZ÷2, NX÷2)
    end
    position += POFFSET
    signal = Signal(signature, position)
end


# starting pressure field
begin
    P0 = zero(v0)
    #P0[1,1] = 1      # endon test
    #P0[1, NX÷2] = 1  # split test
end


# padding arrays, allocating 2 more times for pressure field
begin
    _v = pad_extremes(v0, TAPER)
    _P = pad_zeros_add_zeros_axis(P0, TAPER+1, 3)
end


@time propagate_absorb(grid, _P, _v,  signal)
#@time propagate(grid, _P, _v,  signal, S)


#using BenchmarkTools
## runs 6-7 times, be careful
#print("Propagate absorb ")
#@btime propagate_absorb($grid, $_P, $_v,  $signal)
#print("Propagate pure ")
#@btime propagate_pure($grid, $_P, $_v,  $signal, $S)




using PyPlot

#plt.imshow(_P[1+TAPER:end-TAPER, 1+TAPER:end-TAPER, 1])
#plt.imshow(_P[:, :, 1]; vmin=-.2, vmax=.2)
plt.imshow(_P[:, :, 1])

plt.colorbar()
plt.show()
