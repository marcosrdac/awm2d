using InteractiveUtils
include("./reff.jl")

println("Number of threads = $(nthreads())")


const TAPER = 1


# mesh and solving parameters
begin
    h  = 1.0 # km
    Δt = .001 # s
    NX = 1
    NZ = 1
    NT = 2600
    grid = FDM_Grid(h, Δt, NZ, NX, NT, TAPER)
end


# three layered model parameters
begin
    H1 = (2*NZ)÷8
    H2 = (3*NZ)÷8
    V1 = 3.
    V2 = 4.
    V3 = 5.

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
    signal = rickerwave(ν, Δt)
    array = "center"


    if array === "split"
        S = CartesianIndex(1, NX÷2)
    elseif array === "endon"
        S = CartesianIndex(1, 1)
    elseif array === "center"
        S = CartesianIndex(NZ÷2, NX÷2)
    end

    OFFSET = CartesianIndex(TAPER+∇²r, TAPER+∇²r)
    S += OFFSET
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

@time propagate_absorb(grid, _P, _v,  signal, S)
#@time propagate(grid, _P, _v,  signal, S)


#using BenchmarkTools
## runs 6-7 times, be careful
#print("Propagate absorb ")
#@btime propagate_absorb($grid, $_P, $_v,  $signal, $S)
#print("Propagate pure ")
#@btime propagate($grid, $_P, $_v,  $signal, $S)


#using PyPlot
####imshow(_P[1+TAPER:end-TAPER, 1+TAPER:end-TAPER, 1])
#imshow(_P[:, :, 1])
#plt.show()
