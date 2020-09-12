#!/usr/bin/env julia

using Base.Threads
using BenchmarkTools

macro println() :(println()) end

"including \"propagate_module\", using Propagate module" |> println
include("./propagate_module.jl")
using .Propagate

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

"defining signal" |> println
(sz, sx) = sourceposition(array, NZ, NX)
sourcesignature = discarray(sourcesignaturefile)

shotssignals = [[signal1d(sz, sx, sourcesignature)], [signal1d(sz, sx+20, sourcesignature)]]
# shotssignals = [[signal1d(sz, sx, sourcesignature)], ]

"defining velocity model" |> println
v = discarray(vfile)

"defining initial pressure field" |> println
# P0 = zero(v)

@println

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; filename=$Pfile, stencilorder=2)
@time propagateshots(grid, v, shotssignals, nrec, Δxrec;
                     Pfile=Pfile, seisfile=seisfile, multiseisfile=multiseisfile, stencilorder=8)

run(`python view.py $multiseisfile`)
