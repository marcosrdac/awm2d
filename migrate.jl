#!/usr/bin/env julia

using Base.Threads

macro println() :(println()) end

"including \"propagate_module\", using Propagate module" |> println
@time include("./propagate_module.jl")
@time using .Propagate

@println

"including parameters file" |> println
@time include("./parameters.jl")

@println

"defining grid" |> println
@time grid = FDM_Grid(h, Δt, NZ, NX, NT)

@println

"defining signal" |> println
@time signature = rickerwave(ν, Δt)
position = sourceposition(array, NZ, NX)
signal = Signal1D(signature, position)

@println

"defining velocity model" |> println
# reflectors position
z1 = H1
z2 = z1+H2
# actually defining velocity field
v = Array{Float64}(undef, (NZ, NX))
v[   1:z1,  1:end] .= V1
v[z1+1:z2,  1:end] .= V2
v[z2+1:end, 1:end] .= V3

@println

"defining initial pressure field" |> println
@time P0 = zero(v)

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
