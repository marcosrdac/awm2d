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

# println(1)
# for shots in shotssignals, shot in shots
    # display(shot.x)
    # println()
# end
# println(nx)


#"defining receptors" |> println
#Δxrec = 10
#nrec  = 2*ceil(Int, Δxshot/10)
#shotsrecpositions = [receptorpositions("split", signals[1].x,  Δxrec, nrec)
#                     for signals in shotssignals]
#
## println(1)
## for shots in shotsrecpositions, shot in shots
#    # display(shot)
#    # println()
## end
## println(nx)

shotsrecpositions = [edgeindices(nz, nx) for shot in shotssignals]


"source signal propagation" |> println
# using BenchmarkTools
# @btime propagate($grid, $v, $signal; filename=$Pfile, stencilorder=2)

@time propagateshots(grid, v;
                     seisfile=seisfile,
                     multiseisfile=multiseisfile,
                     shotssignals=shotssignals,
                     shotsrecpositions=shotsrecpositions,
                     )

run(`python view.py $multiseisfile`)
