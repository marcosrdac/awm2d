#!/usr/bin/env julia
using Base.Threads

"including \"propagate_module\", using Propagate module" |> println
include("./propagate_module.jl")
using .Propagate

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDM_Grid(Δz, Δx, Δt, NZ, NX, NT)

"defining velocity model" |> println
v = discarray(v_file)

"defining initial pressure field" |> println
P0 = zero(v)

"defining signal" |> println
source_signature = discarray(source_signature_file)
(sz, sx) = sourceposition(array, NZ, NX)
signal = signal1d(sz, sx, source_signature)

println()

"source signal propagation" |> println
@time propagate_save(grid, P0, v, signal;
                     filename=P_file, stencil_order=stencil_order)


println()


"source signal direct wave propagation" |> println
@time propagate_save_seis(grid, P0, v, signal;
                          filename=direct_seis_file, stencil_order=stencil_order,
                          direct_only=true)


println()


"reverse wave propagation" |> println
begin
    seis = P2seis(discarray(P_file))
    direct_seis = discarray(direct_seis_file)

    # removing direct_wave from seismogram
    @threads for I in eachindex(direct_seis)
        seis[I] -= direct_seis[I]
    end

    seis_signals = seis2signals(seis)
end
@time propagate_save(grid, P0, v, seis_signals; filename=reversed_P_file, stencil_order=stencil_order)


println()


"applying image condition" |> println
@time image_condition(P_file, reversed_P_file, migrated_file)
