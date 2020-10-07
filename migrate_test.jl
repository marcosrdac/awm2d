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

nshot = 9
Δxshot = nx ÷ (1+nshot+4)

shotssignals = [[Signal1D(1, x, signature)]
                for i in (1:nshot) .- (1+nshot÷2)
                for x = nx÷2+i*Δxshot]

"defining receptors" |> println
#Δxrec = 10
#nrec  = 2*ceil(Int, Δxshot/10)
#shotsrecpositions = [receptorpositions("split", signals[1].x,  Δxrec, nrec)
#                     for signals in shotssignals]
shotsrecpositions = [edgeindices(nz, nx) for shot in shotssignals]

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
