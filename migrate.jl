#!/usr/bin/env julia

"including module files" |> println
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
grid = FDMGrid(Δz, Δx, Δt, nz, nx, nt)

"defining signal" |> println
sourcesignature = discarray(sourcesignaturefile)
(sz, sx) = sourceposition(array, nz, nx)
signal = Signal1D(sz, sx, sourcesignature)

"source signal propagation" |> println
@time propagate(grid, v, signal; filename=Pfile)

"source signal direct wave propagation" |> println
@time propagate(grid, v, signal; filename=directseisfile,
                        directonly=true)

"reverse wave propagation" |> println
using Base.Threads
begin
    seis = P2seis(discarray(Pfile))
    directseis = discarray(directseisfile)

    # removing direct_wave from seismogram
    @threads for I in eachindex(directseis)
        seis[I] -= directseis[I]
    end

    seissignals = seis2signals(seis)
end
@time propagate(grid, v, seissignals; filename=reversedPfile)

"applying image condition" |> println
@time imagecondition(Pfile, reversedPfile, migratedfile)
