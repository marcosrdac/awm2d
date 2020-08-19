module Propagate
    using Base.Threads
    using Parameters
    using SparseArrays
    using StaticArrays
    using Mmap: mmap

    # structures
    export FDM_Grid, Signal1D, Signal2D
    # functions
    export gen_3lay_v, rickerwave, sourceposition, discarray, todiscarray, slice_seismogram, image_condition
    export propagate, propagate_save, propagate_save_seis, propagate_2d, propagate_2d_save, propagate_2d_save_seis


    """
        TAPER
    Tapering layer length.
    """
    const TAPER = 60


    """
        ATTENUATION_COEFICIENT
    Positive coefficient value for attenuation in taper.
    """
    const ATTENUATION_COEFICIENT = 0.0035  # previously found
    #const ATTENUATION_COEFICIENT = 0.008    # Átila's
    

    function central_difference_coefs(degree::Integer, order::Integer)
        @assert order % 2 === 0
        p = Int(floor((order+1)/2))
        # defining P matrix
        P = Array{Float64}(undef, 2p+1, 2p+1)
        P[1,:] .= 1
        P[2,:] = -p:p
        for i in 2:size(P, 1)
            P[i,:] = P[2,:] .^ (i-1)
        end
        # defining d matrix
        d = zeros(size(P, 1))
        d[degree+1] = factorial(degree)
        # solving equation P c = d
        c = P\d
    end 

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
    Structure for defining a 2D finite-difference grid.
    """
    struct FDM_Grid{T}
        h::T
        Δt::T
        nz::Int
        nx::Int
        nt::Int
    end


    function gen_3lay_v(nz, nx, h1, h2, V1, V2, V3)
        # three layered model parameters
        h1 = (1*nz)÷3
        h2 = (1*nz)÷3
        # reflectors position
        z1 = h1
        z2 = z1+h2
        # defining velocity field
        v = Array{Float64}(undef, (nz, nx))
        v[   1:z1,  1:end] .= V1
        v[z1+1:z2,  1:end] .= V2
        v[z2+1:end, 1:end] .= V3
        return v
    end


    """
        Signal(signature::AbstractArray{T}, position::CartesianIndex{2})
    Structure for defining a source signal.
    """
    abstract type AbstractSignal end

    struct Signal1D{T} <: AbstractSignal
        signature::AbstractArray{T, 1}
        position::CartesianIndex{2}
    end

    struct Signal2D{T} <: AbstractSignal
        signature::AbstractArray{T, 2}
        position::CartesianIndex{2}
    end

    
    function sourceposition(arrayname::String, NZ::Integer, NX::Integer)
        if arrayname === "split"
            position = CartesianIndex(1, NX÷2+1)
        elseif arrayname === "endon"
            position = CartesianIndex(1, 1)
        elseif arrayname === "center"
            position = CartesianIndex(NZ÷2+1, NX÷2+1)
        end
        position
    end


    function discarray(filename::String, mode::String="r", type::DataType=Float64, dims=())
        if occursin("w", mode)
            @assert dims !== ()
            @assert eltype(dims) <: Int64
            io = open(filename, mode)
            _ndims = length(dims)
            write(io, _ndims, dims...)
            A = mmap(io, Array{type, _ndims}, dims)
            close(io)
        else
            io = open(filename, mode)
            _ndims = read(io, Int64)
            dims = Tuple(read(io, Int64) for i in 1:_ndims)
            A = mmap(io, Array{type, _ndims}, dims)
            close(io)
        end
        return(A)
    end

    function todiscarray(filename::String, A::AbstractArray)
        discA = discarray(filename, "w+", eltype(A), size(A))
        discA .= A
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
        _A = zeros(eltype(A), (size(A, 1)+2*padding,
                            size(A, 2)+2*padding,
                            dim))
        _A[1+padding:end-padding, 1+padding:end-padding, 1] .= A
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
        @assert ν < 0.2 * 1.0/(2.0*Δt)
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
    function h²∇²(A, IL, dz, dx, ∇²_stencil)
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
    function new_p(IP, Iv, cur_P, old_P, v, dz, dx, ∇²_stencil)
        2*cur_P[IP] - old_P[IP] + v[Iv]^2 * Δtoh² * h²∇²(cur_P, IP, dz, dx, ∇²_stencil)
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

    function get_attenuation_factors(TAPER)
        factors = Array{Float64, 1}(undef, TAPER)
        for dist in eachindex(factors)
            factors[dist] = attenuation_factor(dist)
        end
        return(factors)
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


    function offset_signal(signal::AbstractSignal, offset)
        _signal = typeof(signal)(signal.signature,
                                signal.position + offset)
    end


    # defining propagation functions
    macro putsignal(func)
        # adding signal to field before iteration
        funcname = String(func)
        if occursin("2d", funcname)
            quote
                if T <= size(_signal.signature, 1)
                    cur_P[POFFSET+signal.position[1], POFFSET+1:POFFSET+nx] =
                        _signal.signature[1+size(_signal.signature, 1)-T,:]
                end
            end
        else
            quote
                if T <= length(_signal.signature)
                    cur_P[_signal.position] = _signal.signature[T]
                end
            end
        end |> esc
    end

    macro initialize_saved_arrays(func)
        funcname = String(func)
        if occursin("save", funcname)
            if occursin("seis", funcname)
                quote
                    saved_seis = discarray(filename, "w+",Float64, (nt, nx))
                    saved_seis[1,:] .= P[1,:,1]
                end
            else
                quote
                    saved_P = discarray(filename, "w+", Float64, (nz, nx, nt))
                    saved_P[:,:,1] .= P[:,:,1]
                end
            end
        else
            if occursin("seis", funcname)
                quote
                    saved_seis = Array{Float64, 2}(undef, (nt, nx))
                    saved_seis[1,:] .= P[1,:,1]
                end
            else
                quote
                    saved_P = Array{Float64, 3}(undef, (nz, nx, nt))
                    saved_P[:,:,1] .= P[:,:,1]
                end
            end
        end |> esc
    end

    macro savesnapifsave(func)
        funcname = String(func)
        if occursin("save", funcname)
            if occursin("seis", funcname)
                :(saved_seis[T,:] .= P[1,:,new_t])
            else
                :(saved_P[:,:,T] .= P[:,:,new_t])
            end
        end |> esc
    end

    macro returnpropperarrayifnotsave(func)
        funcname = String(func)
        if ! occursin("save", funcname)
            if occursin("seis", funcname)
                :(return(saved_seis))
            else
                :(return(saved_P))
            end
        end |> esc
    end

    for func in (:propagate,    :propagate_save,    :propagate_save_seis,
                 :propagate_2d, :propagate_2d_save, :propagate_2d_save_seis)
        quote
    #==========================================================================#
    function $func(grid, P0, v, signal::AbstractSignal;
                    filename::String="data.bin",
                    stencil_order::Integer=2,
                    direct_only::Bool=false,)
        global TAPER, POFFSET, IPOFFSET
        @unpack h, Δt, nz, nx, nt = grid
        ∇²_stencil = SArray{order+1}(central_difference_coefs(2, order))
        borders = get_taper_sectors(nz, nx, TAPER)
        attenuation_factors = get_attenuation_factors(TAPER)
        Δtoh² = Δt/h^2

        _v = pad_extremes(v, TAPER)
        _signal = offset_signal(signal, IPOFFSET)
        _P = pad_zeros_add_axes(P0, TAPER+1, 3)
        # pressure field of interest
        P = @view _P[1+POFFSET:end-POFFSET, 1+POFFSET:end-POFFSET, :]

        if direct_only
            @views _v[:,TAPER+1:end-TAPER] .= repeat(v[1:1,:], nz+2*TAPER, 1)
        end

        @initialize_saved_arrays($func)  # declares saved_P or saved_seis

        # time loop, order is important
        for T in eachindex(2:nt)
            new_t = mod1(T,   3)
            cur_t = mod1(T-1, 3)
            old_t = mod1(T-2, 3)

            old_P = @view _P[:,:,old_t]
            cur_P = @view _P[:,:,cur_t]
            new_P = @view _P[:,:,new_t]

            # attenuating current and old iteration borders
            for border in borders
                @threads for Iv in border.indices
                    IP = Iv + I∇²r
                    dist = nearest_border_distance(_v, border.id, Iv)
                    cur_P[IP] *= attenuation_factors[dist]
                    old_P[IP] *= attenuation_factors[dist]
                end
            end

            # putting 1D or 2D signal
            @putsignal($func)

            # solve P wave equation
            @threads for Iv in CartesianIndices(_v)
                IP = Iv + I∇²r
                new_P[IP] = new_p(IP, Iv, cur_P, old_P, _v, dz, dx, ∇²_stencil)
            end
            @savesnapifsave($func)
        end
        @returnpropperarrayifnotsave($func)
    end
    #==========================================================================#
        end |> eval
    end

    function slice_seismogram(P)
        seis = copy(P[1,:,:]')
    end

    function image_condition(P_file, reversed_P_file, migrated_file)
        P = discarray(P_file)
        reversed_P = discarray(reversed_P_file)
        (nz, nx, nt) = size(P)
        migrated = discarray(migrated_file, "w+", Float64, (nz, nx))
        migrated .= 0.
        for t in 1:nt
            @views migrated[:,:] .+= P[:,:,t] .* reversed_P[:,:,nt-t+1]
        end
    end
end
