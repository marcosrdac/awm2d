#!/usr/bin/env julia

"including \"propagate_module\", using Propagate module" |> println
include("src/discarrays.jl")
include("src/acoustics2d.jl")
using .Discarrays
using .Acoustics2D

"including parameters file" |> println
include("./parameters.jl")

"defining velocity model" |> println
v = discarray(vfile)
nz, nx = size(v)

"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

"defining sources" |> println
(sz, sx) = sourceposition(array, 1, NX)
sourcesignature = discarray(sourcesignaturefile)

shotssignals = [[Signal1D(sz, sx-20, sourcesignature)],
                [Signal1D(sz, sx,    sourcesignature)],
                [Signal1D(sz, sx+20, sourcesignature)]]

"defining receptors" |> println
shotsrecpositions = [receptorpositions("endon2right", signals[1].x, 5, 10)
                     for signals in shotssignals]

"defining initial pressure field" |> println
# P0 = zero(v)


"source signal propagation" |> println
# using BenchmarkTools
# @btime propagate($grid, $v, $signal; filename=$Pfile, stencilorder=2)

@time propagateshots(grid, v, shotssignals, shotsrecpositions;
                   seisfile=seisfile,
                   multiseisfile=multiseisfile,)

run(`python view.py $multiseisfile`)
