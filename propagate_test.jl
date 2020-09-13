#!/usr/bin/env julia

using Base.Threads
using BenchmarkTools

macro println() :(println()) end

"including \"propagate_module\", using Propagate module" |> println
include("./discarrays.jl")
include("./acoustics2d.jl")
using .Discarrays
using .Acoustics2D

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

"defining signal" |> println
(sz, sx) = sourceposition(array, NX)
sourcesignature = discarray(sourcesignaturefile)
signal = signal1d(sz, sx, sourcesignature)

"defining velocity model" |> println
v = discarray(vfile)

"defining initial pressure field" |> println
P0 = zero(v)
# P0 = collect(1:reduce(*,size(v)))
# P0 = reshape(P0, size(v))

@println

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; Pfile=$Pfile, stencilorder=2)

# @time propagate(grid, v, signal, P0; Pfile=Pfile, stencilorder=16)
# run(`python view.py $Pfile 500`)

# @time propagate(grid, v, signal, P0; seisfile=seisfile, stencilorder=8)

propagate(grid, v, [signal], P0; seisfile=seisfile, stencilorder=2)

run(`python view.py $seisfile 200`)
