#!/usr/bin/env julia

using Base.Threads

macro println() :(println()) end

"including \"propagate_module\", using Propagate module" |> println
include("./propagate_module.jl")
using .Propagate

"including parameters file" |> println
include("./parameters.jl")

"defining grid" |> println
grid = FDM_Grid(h, Δt, NZ, NX, NT)

"defining signal" |> println
position = sourceposition(array, NZ, NX)
source_signature = discarray(source_signature_file)
signal = Signal1D(source_signature, position)

"defining velocity model" |> println
v = discarray(v_file)

"defining initial pressure field" |> println
P0 = zero(v)

@println

"source signal propagation" |> println
@time propagate_save(grid, P0, v, signal; filename=P_file)


@println


"source signal direct wave propagation" |> println
@time propagate_save_seis(grid, P0, v, signal; filename=direct_seis_file, direct_only=true)


@println


"reverse wave propagation" |> println
begin
    seis = slice_seismogram(discarray(P_file))
    direct_seis = discarray(direct_seis_file)

    # removing direct_wave from seismogram
    @threads for I in eachindex(direct_seis)
        seis[I] -= direct_seis[I]
    end

    seis_signal = Signal2D(seis, CartesianIndex(1, 1))
end
@time propagate_2d_save(grid, P0, v, seis_signal; filename=reversed_P_file)


@println


"applying image condition" |> println
@time image_condition(P_file, reversed_P_file, migrated_file)
