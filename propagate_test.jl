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

"defining signal" |> println
(sz, sx) = 1, nx÷2
# (sz, sx) = (100,1)
# (sz, sx) = (250, 178)
sourcesignature = discarray(sourcesignaturefile)
signal = Signal1D(sz, sx, sourcesignature)

"defining initial pressure field" |> println
P0 = zero(v)
# P0 = reshape(1.0:reduce(*,size(v)), size(v))
# P0 = repeat

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; Pfile=$Pfile, stencilorder=2)

# snaps
# @time P = propagate(grid, v, signal, P0; Pfile=Pfile, stencilorder=8)
@time P = propagate_rem(grid, v, signal, P0; Pfile=Pfile, stencilorder=2)
run(`python view.py $Pfile $(300-1)`)

# and seismogram
# @time todiscarray(seisfile, P2seis(P))
# run(`python view.py $seisfile 500`)


# only seis
# @time propagate(grid, v, signal, P0; seisfile=seisfile, stencilorder=8)
# run(`python view.py $seisfile 200`)
