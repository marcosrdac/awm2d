#!/usr/bin/env julia

using Base.Threads
using BenchmarkTools

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
@btime propagate($grid, $P0, $v, $signal; filename=$P_file)
