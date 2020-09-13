#!/usr/bin/env julia

using BenchmarkTools

"including module files" |> println
include("src/discarrays.jl")
include("src/acoustics2d.jl")
using .Discarrays
using .Acoustics2D

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

"defining signal" |> println
(sz, sx) = sourceposition(array, NX)
sourcesignature = discarray(sourcesignaturefile)
signal = Signal1D(sz, sx, sourcesignature)

"defining velocity model" |> println
v = discarray(vfile)

"defining initial pressure field" |> println
P0 = zero(v)
# P0 = reshape(1.0:reduce(*,size(v)), size(v))
# P0 = repeat

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; Pfile=$Pfile, stencilorder=2)

@time propagate(grid, v, signal, P0; Pfile=Pfile, stencilorder=8)
run(`python view.py $Pfile 500`)

# @time propagate(grid, v, signal, P0; seisfile=seisfile, stencilorder=8)

# propagate(grid, v, signal, P0; seisfile=seisfile, stencilorder=8)
# run(`python view.py $seisfile 200`)
