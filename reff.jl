using Base.Threads
using Mmap: mmap
using Parameters
using StaticArrays

"""
    TAPER
The length of the tapering layer used.
"""
const TAPER = 60


"""
    ATTENUATION_COEFICIENT
The positive attenuation coefficient used in exponential formula for getting 
attenuation factors for pressure field modeling.
"""
const ATTENUATION_COEFICIENT = 0.0035


"""
    ∇²
The nine-stencil finite difference laplacian operator over a gaussian curve.
"""
const ∇² = @SArray ([1/6.   4/6.  1/6.
                     4/6. -20/6.  4/6.
                     1/6.   4/6.  1/6.])


"""
    ∇²r
The radius of the laplacian operator used. Useful for padding arrays or getting
correct offsets.
"""
const ∇²r = size(∇², 1)÷2


"""
    I∇²r
The radius of the laplacian operator used, in cartesian index. Useful for 
offsetting positions.
"""
const I∇²r = CartesianIndex(∇²r, ∇²r)


"""
    I1
Cartesian index pointing to [1, 1]. Useful for offsetting positions.
"""
const I1 = CartesianIndex(1, 1)


"""
    POFFSET
    Cartesian index of the [0, 0] position in padded  pressure field. Useful for
offsetting positions.
"""
const POFFSET = TAPER+∇²r


"""
    IPOFFSET
    Cartesian index of the [0, 0] position in padded  pressure field. Useful for
offsetting positions.
"""
const IPOFFSET = CartesianIndex(POFFSET, POFFSET)


"""
    FDM_Grid(h::T, Δt::T, nz::Int, nx::Int, nt::Int)
Structure for defining a 2D finite differences grid.
"""
struct FDM_Grid{T}
    h::T
    Δt::T
    nz::Int
    nx::Int
    nt::Int
end


"""
    Signal(signature::AbstractArray{T}, position::CartesianIndex{2})
Structure for defining a source signal.
"""
struct Signal{T}
    signature::AbstractArray{T}
    position::CartesianIndex{2}
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
function pad_zeros_add_axes(A::AbstractArray,
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


"""
    h²∇²(A, IL)
Function get the laplacian of a point inside an array. In order to optimize code,
it actually calculates the laplacian without actually dividing the values by h².
Then it's a h² laplacian function.
"""
function h²∇²(A, IL)
    global ∇², ∇²r, I∇²r, I1
    result = zero(eltype(A))
    @inbounds for I in CartesianIndices(∇²)
        J = IL - I∇²r - I1 + I
        result += A[J] * ∇²[I]
    end
    return result
end


"""
    new_p(IP, Iv, cur_P, old_P, v, Δtoh²)
Function that calculates the next P value at point IP in pressure field,
considering propagation velocity at that point v[Iv]; cur_P and old_P are,
respectively, the current pressure field and old pressure field arrays; Δtoh²
is the relation Δt/h² used in the modeling process.
"""
function new_p(IP, Iv, cur_P, old_P, v, Δtoh²)
    2*cur_P[IP] - old_P[IP] + v[Iv]^2 * Δtoh² * h²∇²(cur_P, IP)
end


"""
    function nearest_border_distance(A, border, R)
Add description...
"""
function nearest_border_distance(A, border, R)
    # orthogonal
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


"""
    attenuation_factor(dist)
Add description...
"""
function attenuation_factor(dist)
    global ATTENUATION_COEFICIENT, TAPER
    # for testing in python return(exp(-(0.0035*x))
    exp(-(ATTENUATION_COEFICIENT*(TAPER-dist))^2)
end


"""
    propagate_pure(grid, P, v, signal)
Add description...
"""
function propagate_pure(grid, P, v, signal)
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


"""
    Sector(id::T, indices::CartesianIndices{2, Tuple{UnitRange{T},
                                                    UnitRange{T}}})
Add description...
"""
struct Sector{T}
    id::T
    indices::CartesianIndices{2, Tuple{UnitRange{T}, UnitRange{T}}}
end


"""
    get_tapered_sectors(nz, nx, taper=0)
Add description...
"""
function get_tapered_sectors(nz, nx, taper=0)
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

function get_taper_sectors(nz, nx, taper)
    return get_tapered_sectors(nz, nx, taper)[(1:end .!= 5)]
end


function offset_signal(signal::Signal, offset)
    _signal = Signal(signal.signature,
                     signal.position + offset)
end

"""
    propagate(grid, P, v, signal)
Add description...
"""
function propagate(grid, P0, v, signal)
    global TAPER, POFFSET, IPOFFSET
    @unpack h, Δt, nz, nx, nt = grid
    Δtoh² = Δt/h^2

    _v = pad_extremes(v, TAPER)
    _signal = offset_signal(signal, IPOFFSET)
    _P = pad_zeros_add_axes(P0, TAPER+1, 3)

    borders = get_taper_sectors(nz, nx, TAPER)

    # time loop, order is important
    for timeiter in eachindex(1:nt)
        new_t = mod1(timeiter+1, 3)
        cur_t = mod1(timeiter,   3)
        old_t = mod1(timeiter-1, 3)

        old_P = @view _P[:,:,old_t]
        cur_P = @view _P[:,:,cur_t]
        new_P = @view _P[:,:,new_t]

        # attenuating current and old iteration borders
        for border in borders
            @threads for Iv in border.indices
                IP = Iv + I∇²r
                dist = nearest_border_distance(_v, border.id, Iv)
                att = attenuation_factor(dist)
                cur_P[IP] *= att
                old_P[IP] *= att
            end
        end

        # adding signal to field before iteration
        if timeiter <= length(_signal.signature)
            cur_P[_signal.position] = _signal.signature[timeiter]
        end


        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(_v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(IP, Iv, cur_P, old_P, _v, Δtoh²)
        end
    end

    end_t = mod1(nt+1, 3)
    return(_P[1+(TAPER+1):end-(TAPER+1), 1+(TAPER+1):end-(TAPER+1), end_t])
    #return(_P[:,:,end_t])
end


function propagate_save(grid, P0, v, signal;
                        filename::String,
                        only_seis::Bool=false,
                        return_value=false)
    global TAPER, POFFSET, IPOFFSET
    @unpack h, Δt, nz, nx, nt = grid
    borders = get_taper_sectors(nz, nx, TAPER)
    Δtoh² = Δt/h^2

    _v = pad_extremes(v, TAPER)
    _signal = offset_signal(signal, IPOFFSET)
    _P = pad_zeros_add_axes(P0, TAPER+1, 3)
    # pressure field of interest
    P = @view _P[1+POFFSET:end-POFFSET, 1+POFFSET:end-POFFSET, :]

    if only_seis
        saved_dims = (nt, nx)
    else
        saved_dims = (nz, nx, nt)
    end
    saved_ndims = length(saved_dims)

    # setting up disk array for saving output (P_saved)
    io = open(filename, "w+")
    write(io, saved_ndims, saved_dims...)

    if only_seis
        saved_seis = mmap(io, Array{Float64, saved_ndims}, saved_dims)
        saved_seis[1,:] .= P[1,:,1]
    else
        saved_P = mmap(io, Array{Float64, saved_ndims}, saved_dims)
        saved_P[:,:,1] .= P[:,:,1]
    end

    close(io)

    # time loop, order is important
    for timeiter in eachindex(2:nt)
        new_t = mod1(timeiter,   3)
        cur_t = mod1(timeiter-1, 3)
        old_t = mod1(timeiter-2, 3)

        old_P = @view _P[:,:,old_t]
        cur_P = @view _P[:,:,cur_t]
        new_P = @view _P[:,:,new_t]

        # attenuating current and old iteration borders
        for border in borders
            @threads for Iv in border.indices
                IP = Iv + I∇²r
                dist = nearest_border_distance(_v, border.id, Iv)
                att = attenuation_factor(dist)
                cur_P[IP] *= att
                old_P[IP] *= att
            end
        end

        # adding signal to field before iteration
        if timeiter <= length(_signal.signature)
            cur_P[_signal.position] = _signal.signature[timeiter]
        end

        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(_v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(IP, Iv, cur_P, old_P, _v, Δtoh²)
        end

        if only_seis
            saved_seis[timeiter,:] .= P[1,:,new_t]
        else
            saved_P[:,:,timeiter] .= P[:,:,new_t]
        end
    end
end


"""
    save_seis(grid, P, v, signal)
Add description...
"""
function save_seis(grid, P0, v, signal; filename::String)
    global TAPER, POFFSET, IPOFFSET
    @unpack h, Δt, nz, nx, nt = grid
    Δtoh² = Δt/h^2
    borders = get_taper_sectors(nz, nx, TAPER)
    seis = Array{Float64, 2}(undef, (nt, nx))

    _v = pad_extremes(v, TAPER)
    _signal = offset_signal(signal, IPOFFSET)
    _P = pad_zeros_add_axes(P0, TAPER+1, 3)

    # differs
    io = open(filename, "w")

    write(io, size(seis)...)

    cur_P = @view _P[:,:, 1]
    seis[1, :] = cur_P[1+POFFSET, 1+POFFSET:end-POFFSET]

    # time loop, order is important
    for timeiter in eachindex(1:nt)
        new_t = mod1(timeiter+1, 3)
        cur_t = mod1(timeiter,   3)
        old_t = mod1(timeiter-1, 3)

        old_P = @view _P[:,:,old_t]
        cur_P = @view _P[:,:,cur_t]
        new_P = @view _P[:,:,new_t]

        # attenuating current and old iteration borders
        for border in borders
            @threads for Iv in border.indices
                IP = Iv + I∇²r
                dist = nearest_border_distance(_v, border.id, Iv)
                att = attenuation_factor(dist)
                cur_P[IP] *= att
                old_P[IP] *= att
            end
        end

        # adding signal to field before iteration
        if timeiter <= length(_signal.signature)
            cur_P[_signal.position] = _signal.signature[timeiter]
        end


        # spacial loop, order is not important
        @threads for Iv in CartesianIndices(_v)
            IP = Iv + I∇²r
            new_P[IP] = new_p(IP, Iv, cur_P, old_P, _v, Δtoh²)
        end

        # differs
        seis[timeiter, :] = new_P[1+POFFSET,
                                  1+POFFSET:end-POFFSET]
    end

    write(io, seis)
    close(io)
end
