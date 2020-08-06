#!/usr/bin/env julia

include("./parameters.jl")
using Base.Threads


"source signal propagation" |> println
@time propagate_save(grid, P0, v, signal; filename=P_file)


println()


"source signal direct wave propagation" |> println
@time propagate_save_seis(grid, P0, v, signal; filename=direct_seis_file, direct_only=true)


println()


"reverse wave propagation" |> println
begin
    seis = slice_seismogram(hddarray(P_file))
    direct_seis = hddarray(direct_seis_file)

    # removing direct_wave from seismogram
    @threads for I in eachindex(direct_seis)
        seis[I] -= direct_seis[I]
    end

    seis_signal = Signal2D(seis, CartesianIndex(1, 1))
end
@time propagate_2d_save(grid, P0, v, seis_signal; filename=reversed_P_file, direct_only=true)


println()


"applying image condition" |> println
@time image_condition(P_file, reversed_P_file, migrated_file)
