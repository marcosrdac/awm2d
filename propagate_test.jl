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
# uncomment below to ignore parameter file's source position
# (sz, sx) = (1, 1)  # (nz÷2, nx÷2)
sourcesignature = discarray(sourcesignaturefile)
signal = Signal1D(sz, sx, sourcesignature)

"defining initial pressure field" |> println
P0 = zero(v)

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; Pfile=$Pfile, stencilorder=2)


# === Modeling snaps === #
@time P = propagate(grid, v, signal, P0; rem=true, pseudo=false, Pfile=Pfile, stencilorder=8)
# run(`python view.py $Pfile $(nt-1)`)
# --- also save only seismogram to seisfile --- #
# @time todiscarray(seisfile, P2seis(P))
# run(`python view.py $seisfile 500`)


# === Modeling only seismogram === #
# Modeling only seismogram
# @time propagate(grid, v, signal, P0; seisfile=seisfile, stencilorder=8)
# run(`python view.py $seisfile 200`)
