#!/usr/bin/env julia

"including \"propagate_module\", using Propagate module" |> println
include("src/discarrays.jl")
include("src/acoustics2d.jl")
using .Discarrays
using .Acoustics2D

# "including parameters file" |> println
# include("./parameters.jl")

# "defining grid" |> println
# grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

# "defining signal" |> println
# (sz, sx) = sourceposition(array, 1, NX)
# sourcesignature = discarray(sourcesignaturefile)

# shotssignals = [[signal1d(sz, sx-20, sourcesignature)], [signal1d(sz, sx, sourcesignature)], [signal1d(sz, sx+20, sourcesignature)]]
# # shotssignals = [[signal1d(sz, sx, sourcesignature)], ]

# "defining velocity model" |> println
# v = discarray(vfile)

# "defining initial pressure field" |> println
# # P0 = zero(v)

# @println

defini

#"source signal propagation" |> println
# using BenchmarkTools
## @btime propagate($grid, $v, $signal; filename=$Pfile, stencilorder=2)
#@time propagateshots(grid, v, shotssignals, nrec, Δxrec;
#                     Pfile=Pfile, seisfile=seisfile, multiseisfile=multiseisfile, stencilorder=8)
#
#run(`python view.py $multiseisfile`)
