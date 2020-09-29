#!/usr/bin/env julia

"including module files" |> println
include("src/discarrays.jl")
include("src/acoustics2d.jl")
using .Discarrays
using .Acoustics2D
"including parameters file" |> println
include("parameters.jl")


"defining velocity model" |> println
v = discarray(vfile)
nz, nx = size(v)
"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, nz, nx, nt)


"defining sources" |> println
signature = discarray(sourcesignaturefile)

(sz1, sx1) = 1, nx÷2 - (nx÷4+1)
(sz2, sx2) = 1, nx÷2 - (nx÷8+1)
(sz3, sx3) = 1, nx÷2
(sz4, sx4) = 1, nx÷2 + (nx÷8+1)
(sz5, sx5) = 1, nx÷2 + (nx÷4+1)

shotssignals = [[Signal1D(sz1, sx1, signature)],
                [Signal1D(sz2, sx2, signature)],
                [Signal1D(sz3, sx3, signature)],
                [Signal1D(sz4, sx4, signature)],
                [Signal1D(sz5, sx5, signature)]]

"defining receptors" |> println
#                                             array             sx  Δx  nrec
shotsrecpositions = [receptorpositions("split", signals[1].x,  ceil(Int, 5/321*nx),  12)
                     for signals in shotssignals]


"source signal propagation" |> println
# using BenchmarkTools
# @btime propagate($grid, $v, $signal; filename=$Pfile, stencilorder=2)

@time migrate(
        grid, v;
        shotssignals=shotssignals,
        Pfile=Pfile,
        multiseisfile=multiseisfile,
        reversedPfile=reversedPfile,
        migratedfile=migratedfile,
        shotsrecpositions=shotsrecpositions,
)

run(`python view.py $migratedfile`)
