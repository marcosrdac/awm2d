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
grid = FDM_Grid(Δz, Δx, Δt, NZ, NX, NT)

"defining signal" |> println
(sz, sx) = sourceposition(array, NX)
source_signature = discarray(source_signature_file)
signal = signal1d(sz, sx, source_signature)

"defining velocity model" |> println
v = discarray(v_file)

"defining initial pressure field" |> println
P0 = zero(v)

@println

"source signal propagation" |> println
# @btime propagate($grid, $v, $signal; P_file=$P_file, stencil_order=2)
@time propagate(grid, v, signal; Pfile=P_file, stencil_order=2)
# @time propagate(grid, v, signal; onlyseis=false, stencil_order=2)
