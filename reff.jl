using Base.Threads
using Parameters
using StaticArrays

const ∇² = @SArray ([1/6.   4/6.  1/6.
                     4/6. -20/6.  4/6.
                     1/6.   4/6.  1/6.])

const ∇²r = size(∇², 1)÷2
const I∇²r = CartesianIndex(∇²r, ∇²r)
const I1 = CartesianIndex(1, 1)
const ATTENUATION_COEFICIENT = 0.0035


struct FDM_Grid{T}
    h::T
    Δt::T
    nz::Int
    nx::Int
    nt::Int
    taper::Int
end

struct Signal{T}
    signature::AbstractArray{T}
    position::CartesianIndex{2}
end

function h²∇²(cur_P, IP)
    global ∇², ∇²r, I∇²r, I1
    result = zero(eltype(cur_P))
    @inbounds for I in CartesianIndices(∇²)
        J = IP - I∇²r - I1 + I
        result += cur_P[J] * ∇²[I]
    end
    return result
end


function new_p(IP, Iv, cur_P, old_P, v, Δtoh²)
    2*cur_P[IP] - old_P[IP] + v[Iv]^2 * Δtoh² * h²∇²(cur_P, IP)
end

function nearest_border_distance(A, border, R)
    if border === 4      R[2]
    elseif border === 8  R[1]
    elseif border === 2  size(A, 1)-R[1]+1
    elseif border === 6  size(A, 2)-R[2]+1
    # diagonal
    elseif border === 7  min(R[1], R[2])
    elseif border === 9  min(size(A, 2)-R[2]+1, R[1])
    elseif border === 1  min(size(A, 1)-R[1]+1, R[2])
    elseif border === 3  min(size(A, 1)-R[1], size(A, 2)-R[2]) + 1
    else zero(eltype(R))
    end
end


# for testing in python return(exp(-(0.0035*x))
function attenuation_factor(taper, dist)
    exp(-(ATTENUATION_COEFICIENT*(taper-dist))^2)
end

# saved of this plague
#function new_attenuated_p(IP, Iv, cur_P, old_P, v, Δtoh², taper, att_coef, border::Integer=1)
#    dist = nearest_border_distance(v, border, Iv)
#    att = attenuation_factor(taper, dist, att_coef)
#
#    att * (2*cur_P[IP] - att * old_P[IP] + v[Iv]^2 * Δtoh² * h²∇²(cur_P, IP))
#end


function propagate(grid, P, v, signal)
    @unpack h, Δt, nz, nx, nt = grid
    Δtoh² = Δt/h^2
    # time loop, order is important
    for timeiter in eachindex(1:nt)
        # which of the thre times of P are what?
        new_t = mod1(timeiter+1, 3)    # [2] now
        cur_t = mod1(timeiter,   3)    # [1] last
        old_t = mod1(timeiter+2, 3)    # [3] before last

        old_P = @view P[:,:,old_t]
        cur_P = @view P[:,:,cur_t]
        new_P = @view P[:,:,new_t]

        # adding signal to field before iteration
        if timeiter <= size(signal.signature, 1)
            cur_P[signal.position] = signal.signature[timeiter]
        end

        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(IP, Iv, cur_P, old_P, v, Δtoh²)
        end
    end
end


struct Sector{T}
    id::T  # THIS WHAS THE PROBLEM: SETTING id TO Integer
    indices::CartesianIndices{2, Tuple{UnitRange{T}, UnitRange{T}}}
end


function get_sectors(nz, nx, taper=0)
    sectors = Array{Sector}(undef, 9)
    sectors[5] = Sector(5, CartesianIndices((taper+1:taper+nz,
                                             taper+1:taper+nx)))
    sectors[7] = Sector(7, CartesianIndices((1:taper,
                                             1:taper)))
    sectors[4] = Sector(4, CartesianIndices((taper+1:taper+nz,
                                             1:taper)))
    sectors[1] = Sector(1, CartesianIndices((taper+nz+1:nz+2taper,
                                             1:taper)))
    sectors[8] = Sector(8, CartesianIndices((1:taper,
                                             taper+1:taper+nx)))
    sectors[2] = Sector(2, CartesianIndices((taper+nz+1:nz+2taper,
                                             taper+1:taper+nx)))
    sectors[9] = Sector(9, CartesianIndices((1:taper,
                                             taper+nx+1:nx+2taper)))
    sectors[6] = Sector(6, CartesianIndices((taper+1:taper+nz,
                                             taper+nx+1:nx+2taper)))
    sectors[3] = Sector(3, CartesianIndices((taper+nz+1:nz+2taper,
                                             taper+nx+1:nx+2taper)))
    return(sectors)
end



function propagate_absorb(grid, P, v, signal)  # rewrite as source struct
    @unpack h, Δt, nz, nx, nt, taper = grid
    Δtoh² = Δt/h^2
    borders = get_sectors(nz, nx, taper)[(1:end .!= 5)]

    # time loop, order is important
    for timeiter in eachindex(1:nt)
        new_t = mod1(timeiter+1, 3)    # [2] now
        cur_t = mod1(timeiter,   3)    # [1] last
        old_t = mod1(timeiter+2, 3)    # [3] before last

        old_P = @view P[:,:,old_t]
        cur_P = @view P[:,:,cur_t]
        new_P = @view P[:,:,new_t]

        # adding signal to field before iteration
        if timeiter <= size(signal.signature, 1)
            cur_P[signal.position] = signal.signature[timeiter]
        end

        # attenuating current and old iteration borders
        for border in borders
            @threads for Iv in border.indices
                IP = Iv + I∇²r
                dist = nearest_border_distance(v, border.id, Iv)
                att = attenuation_factor(taper, dist)
                cur_P[IP] *= att
                old_P[IP] *= att
            end
        end

        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(IP, Iv, cur_P, old_P, v, Δtoh²)
        end
    end
end





"""
    pad_extremes(A::AbstractArray, padding::Integer)
Pads extremes of A by a determined padding length.
"""
function pad_extremes(A::AbstractArray, padding::Integer)
    _A = Array{eltype(A)}(undef, (size(A, 1)+2*padding,
                             size(A, 2)+2*padding))
    # center
    _A[1+padding:size(A,1)+padding, 1+padding:size(A,2)+padding] .= A
    # borders
    _A[1:padding, 1+padding:end-padding] .= repeat(A[1:1,1:end], padding,1)
    _A[1+padding:end-padding, 1:padding] .= repeat(A[1:end,1:1], 1, padding)
    _A[1+end-padding:end, 1+padding:end-padding] .= repeat(A[end:end,1:end], padding,1)
    _A[1+padding:end-padding, 1+end-padding:end] .= repeat(A[1:end,end:end], 1, padding)
    # diagonal edges
    _A[1:padding, 1:padding] .= A[1,1]
    _A[1+end-padding:end, 1:padding] .= A[end,1]
    _A[1+end-padding:end, 1+end-padding:end] .= A[end,end]
    _A[1:padding, 1+end-padding:end] .= A[1,end]
    _A
end


"""
    pad_extremes(A::AbstractArray, padding::Integer=0, dim::Integer=1)
Pads extremes of A by a determined padding length, adds a dimension. All new elements are initialized as zeros.
"""
function pad_zeros_add_zeros_axis(A::AbstractArray,
                                  padding::Integer=1,
                                  dim::Integer=1)
    _A = zeros(eltype(A), (size(P0, 1)+2*padding,
                           size(P0, 2)+2*padding,
                           dim))
    _A[1+padding:end-padding,
       1+padding:end-padding,
       1] .= A
    _A
end


"""
    norm_ricker(t::Real, σ::Real=1)
Compute the negative normalized ricker function of t with standard deviation σ.
"""
function norm_ricker(t::Real, σ::Real=1)
    2/(√(3*σ)π^(1/4)) * (1-(t/σ)^2) * exp(-t^2/(2σ^2))
end


"""
    rickerwave(ν::Real, Δt::Real)
Compute a ricker wave with main frequency ν, sampled in time intervals Δt.

For a stable signal assert ``ν << \\frac{1}{2 Δt}``.
"""
function rickerwave(ν::Real, Δt::Real)
    @assert ν < 0.02 * 1.0/(2.0*Δt)
    function ricker(t::Real)
        return((1-2*t^2) * exp(-t^2))
    end
    len = 2 * (2.2/(ν*Δt)) ÷ 2
    t = (π*ν*Δt) .* (-len÷2:len÷2)
    ricker.(t)
end


"""
   reduce_pointwise_product(A::AbstractArray, B::AbstractArray)
Performs the pointwise multiplication sum of A over B.
"""
function reduce_pointwise_product(A::AbstractArray, B::AbstractArray)
    prod_sum = zero(eltype(A))
    @fastmath @inbounds @simd for i in eachindex(A)
        prod_sum += A[i] .* B[i]
    end
    prod_sum
end
