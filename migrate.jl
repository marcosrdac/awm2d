#!/usr/bin/env julia
using Base.Threads

"including \"propagate_module\", using Propagate module" |> println
include("./propagate_module.jl")
using .Propagate

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDMGrid(Δz, Δx, Δt, NZ, NX, NT)

"defining velocity model" |> println
v = discarray(vfile)

"defining initial pressure field" |> println
P0 = zero(v)

"defining signal" |> println
sourcesignature = discarray(sourcesignaturefile)
(sz, sx) = sourceposition(array, NZ, NX)
signal = signal1d(sz, sx, sourcesignature)

println()

"source signal propagation" |> println
@time propagatesave(grid, P0, v, signal;
                     filename=Pfile, stencilorder=stencilorder)


println()


"source signal direct wave propagation" |> println
@time propagatesaveseis(grid, P0, v, signal;
                          filename=directseisfile, stencilorder=stencilorder,
                          directonly=true)


println()


"reverse wave propagation" |> println
begin
    seis = P2seis(discarray(Pfile))
    directseis = discarray(directseisfile)

    # removing direct_wave from seismogram
    @threads for I in eachindex(directseis)
        seis[I] -= directseis[I]
    end

    seissignals = seis2signals(seis)
end
@time propagatesave(grid, P0, v, seissignals; filename=reversedPfile, stencilorder=stencilorder)


println()


"applying image condition" |> println
@time imagecondition(Pfile, reversedPfile, migratedfile)
