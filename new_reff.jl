using Base.Threads
using Parameters
using StaticArrays

const ∇² = @SArray ([1/6.   4/6.  1/6.
                     4/6. -20/6.  4/6.
                     1/6.   4/6.  1/6.])

const ∇²r = size(∇², 1)÷2
const I∇²r = CartesianIndex(∇²r, ∇²r)
const I1 = CartesianIndex(1, 1)


struct FDM_Grid{T}
    h::T
    Δt::T
    nz::Int
    nx::Int
    nt::Int
    taper::Int
end


function new_p(cur_P, old_P, v, Δtoh², IP, Iv)
    2*cur_P[IP] - old_P[IP] + v[Iv]^2 * Δtoh² * h²∇²(cur_P, IP)
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


function propagate(grid, P, v, signal, signal_position)
    @unpack h, Δt, nz, nx, nt, taper = grid
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
        if timeiter <= size(signal,1)
            cur_P[signal_position[1], signal_position[2]] = signal[timeiter]
        end

        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(cur_P, old_P, v, Δtoh², IP, Iv)
        end
    end
end





#function new_p(cur_P, old_P, v, Δtoh², zP, xP, zv, xv)
#    2*cur_P[zP, xP] - old_P[zP, xP] + v[zv, xv]^2 * Δt∇²(cur_P, Δtoh², zP, xP)
#end

function nearest_border_distance(nz, nx, border, z, x)
    if border === 4      x
    elseif border === 8  z
    elseif border === 2  nz-z+1
    elseif border === 6  nx-x+1
    # diagonal
    elseif border === 7  min(x, z)
    elseif border === 9  min(nx-x+1, z)
    elseif border === 1  min(x, nz-z+1)
    elseif border === 3  min(nx-x, nz-z)+1
    else                 0
    end
end



function attenuation(taper::Integer, dist::Integer, coef::Real=0.008)
    exp(-(coef*(taper-dist))^2)
end


function new_p_taper(cur_P, old_P, v, Δtoh², zP, xP, zv, xv, border, taper, coef=0.008)
    dist = nearest_border_distance(size(v, 1), size(v, 2), border, zv, xv)
    att = attenuation(taper, dist, coef)
    att * (2*cur_P[zP, xP] -  old_P[zP, xP] + v[zv, xv]^2 * Δt∇²(cur_P, Δtoh², zP, xP))
end



function propagate_etik(P::Array{Float64,3},
                        v::Array{Float64,2},
                        h::Float64,
                        Δt::Float64,
                        nt::Int64,
                        signal::Array{Float64,1}=[],)
    Δtoh² = Δt/h^2

    @inbounds @fastmath for timeiter in eachindex(1:nt)
        new_t = 1 + (timeiter)%3    # 2
        cur_t = 1 + (timeiter+2)%3  # 1
        old_t = 1 + (timeiter+1)%3  # 3

        old_P = @view P[:,:,old_t]
        cur_P = @view P[:,:,cur_t]
        new_P = @view P[:,:,new_t]

        if timeiter <= size(signal,1)
            cur_P[signal_pos[1], signal_pos[2]] = signal[timeiter]
        end

        @inbounds @fastmath @threads for spciter in CartesianIndices(v)
            begin
                zv, xv = spciter[1], spciter[2]
                zP, xP = zv+∇²r, xv+∇²r
                new_P[zP, xP] = new_p(cur_P, old_P, v, Δtoh²,
                                      zP, xP, zv, xv)
            end
        end
    end
end


function propagate_absorb(P::Array{Float64,3},
                          v::Array{Float64,2},
                          h::Float64, Δt::Float64, nt::Int64,
                          signal::Array{Float64,1}=[],
                          taper::Int64=0)
    Δtoh² = Δt/h^2
    nz = size(P, 1)-2*(taper+1)
    nx = size(P, 2)-2*(taper+1)
    if taper === 0
        borders = [5]
        limits = [(taper+1:taper+nz,     taper+1:taper+nx)]
    else
        borders = [5,  7, 4, 1,  8,  2,  9, 6, 3]
        limits = [# center
                     (taper+1:taper+nz,     taper+1:taper+nx),
                  # borders
                           (1:taper,              1:taper),       # top left
                     (taper+1:taper+nz,           1:taper),       # left
                  (taper+nz+1:nz+2taper,          1:taper),       # bottom left
                           (1:taper,        taper+1:taper+nx),   # top
                  (taper+nz+1:nz+2taper,    taper+1:taper+nx),   # bottom
                           (1:taper,     taper+nx+1:nx+2taper),   # top right
                     (taper+1:taper+nz,  taper+nx+1:nx+2taper),   # right
                  (taper+nz+1:nz+2taper, taper+nx+1:nx+2taper),]  # bottom right
    end
              
    @inbounds @fastmath for timeiter in eachindex(1:nt)
        new_t = 1 + (timeiter)%3    # 2
        cur_t = 1 + (timeiter+2)%3  # 1
        old_t = 1 + (timeiter+1)%3  # 3

        old_P = @view P[:,:,old_t]
        cur_P = @view P[:,:,cur_t]
        new_P = @view P[:,:,new_t]

        if timeiter <= size(signal,1)
            cur_P[signal_pos[1], signal_pos[2]] = signal[timeiter]
        end
        #cur_P[signal_pos[1],signal_pos[2]] = (timeiter < size(signal,1) ? signal[timeiter] : 0)

        @inbounds @fastmath @threads for spciter in CartesianIndices(limits[1])
            begin
                zv, xv = spciter[1], spciter[2]
                zP, xP = zv+1, xv+1
                new_P[zP, xP] = new_p(cur_P, old_P, v, Δtoh²,
                                      zP, xP, zv, xv)
            end
        end

        for border in 2:size(borders, 1)
            @inbounds @fastmath @threads for spciter in
                                                CartesianIndices(limits[border])
                begin
                    zv, xv = spciter[1], spciter[2]
                    zP, xP = zv+1, xv+1
                    # new_P[zP, xP] = new_p(cur_P, old_P, v, Δtoh²,
                                        # zP, xP, zv, xv)
                    new_P[zP, xP] = new_p_taper(cur_P, old_P, v, Δtoh²,
                                                zP, xP, zv, xv,
                                                borders[border], taper, 0.008)
                end
            end
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
