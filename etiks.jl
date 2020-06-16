const ∇²r = 1
const h²∇² = fill(1.0, 3, 3)


ABS=0
nz=400
nx=400
nt=4
N = 400
P = rand(N, N, 3)
v = rand(N, N)
signal = rand(N)
nz = N-2
nx = N-2
nt = nt
h = 1.0
Δt = 1.0
ABS = 0
signal_position = [1,1] .+ ABS

############################################################################################
# Idiomatic Julia
############################################################################################
using Base.Threads
using BenchmarkTools
using StaticArrays
using Parameters
using LinearAlgebra

struct GridDefinition{T}
    nz::Int
    nx::Int
    nt::Int
    h::T
    Δt::T
    ABS::Int
    ∇²r::Int
    h²∇²::SArray{Tuple{3, 3}, T, 2, 9}
end

function propagate_idiomatic(P, v, signal, grid_def)
    @unpack nz, nx, ABS, Δt, h, nt = grid_def
    Δtoh² = Δt/h^2
    Plimits = (2:nz+2*ABS+1, 2:nx+2*ABS+1)
    signal_position = SVector(1, 1) .+ ABS

    # time loop, order is important
    for timeiter in eachindex(1:nt)
        # which of the thre times of P are what?
        #new_t =     # [2] now
        #cur_t =     # [1] last
        #old_t =     # [3] before last

        old_P = @view P[:,:,mod1(timeiter+1, 3)]
        cur_P = @view P[:,:,mod1(timeiter,   3)]
        new_P = @view P[:,:,mod1(timeiter+2, 3)]

        # adding signal to field before iteration
        cur_P[signal_position[1],signal_position[2]] =
             (timeiter < size(signal,1) ? signal[timeiter] : zero(eltype(signal)))

        # spacial loop, order is not important
        @threads for I in CartesianIndices(Plimits)
            new_P[I] = new_p_idiomatic(cur_P, old_P, v, Δtoh², I, grid_def)
        end
    end
end

function new_p_idiomatic(cur_P, old_P, v, Δtoh², I, grid_def)
    @unpack h²∇², ∇²r = grid_def  # (constant) 3x3 matrix
    IOFFSET = CartesianIndex(size(h²∇²) .÷ 2)
    h²∇²cur_P = ∇²_kernel(cur_P, h²∇², I, IOFFSET)
    new_p = 2*cur_P[I] - old_P[I] + v[I - CartesianIndex(1, 1)]^2*Δtoh² * h²∇²cur_P
end

function ∇²_kernel(cur_P, h²∇², ZP, IOFFSET)
    result = zero(eltype(cur_P))
    @inbounds for I in CartesianIndices(h²∇²)
        OFFSET = ZP + I - IOFFSET
        cur_result = cur_P[OFFSET] * h²∇²[I]
        result += cur_result
    end
    return result
end

h²∇²_static = SMatrix{3, 3}(fill(1.0, 9))
grid_def = GridDefinition(nz, nx, nt, h, Δt, ABS, ∇²r, h²∇²_static)

print("começando")
@btime propagate_idiomatic($P, $v, $signal, $grid_def)
print("finalizou")
