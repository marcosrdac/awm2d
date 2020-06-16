using BenchmarkTools
using InteractiveUtils
include("./new_reff.jl")

println("Number of threads = $(nthreads())")


const TAPER = 80


# mesh and solving parameters
begin
    h  = 1.0 # km
    Δt = .001 # s
    NX = 250
    NZ = 250
    NT = 432
end


# three layered model parameters
begin
    H1 = (2*NZ)÷6
    H2 = (3*NZ)÷6
    V1 = 3.
    V2 = 3.
    V3 = 3.

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
        signal_pos = [1, NX÷2] .+ (TAPER+∇²r)
    elseif array === "endon"
        signal_pos = [1, 1] .+ (TAPER+1)
    elseif array === "center"
        signal_pos = [NZ÷2, NX÷2] .+ (TAPER+∇²r)    
    end
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


grid = FDM_Grid(h, Δt, NZ, NX, NT, TAPER)
println("foi o grid")
#propagate(grid, _P, _v, signal)
@btime propagate($grid, $_P, $_v,  $signal, $signal_pos)
#@time propagate(grid, _P, _v,  signal)
#println("foi o calc 1")
#@btime propagate($grid, $_P, $_v,  $signal)
#println("foi o calc 2")
##@time propagate(_P, _v, h, Δt, NT, signal)


###using Profile
using PyPlot
####imshow(_P[1+taper:end-taper, 1+taper:end-taper, 1])
imshow(_P[:, :, 1])
plt.show()
#
##using Plots
##heatmap(_P[:, :, 1])
##savefig("test.png")
