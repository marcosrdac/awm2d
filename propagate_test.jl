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
position = sourceposition(array, NZ, NX)
source_signature = discarray(source_signature_file)
signal = Signal1D(source_signature, position)

"defining velocity model" |> println
H1 = (1*NZ)÷3
H2 = (1*NZ)÷3
V1 = 3. # km/s
V2 = 5. # km/s
V3 = 9. # km/s
v = todiscarray(v_file, gen_3lay_v(NZ, NX, H1, H2, V1, V2, V3))

"defining initial pressure field" |> println
P0 = zero(v)

@println

"source signal propagation" |> println
# @btime propagate_save($grid, $P0, $v, $signal; filename=$P_file)
@time propagate(grid, P0, v, signal; filename=P_file)
