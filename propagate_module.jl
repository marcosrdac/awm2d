module Propagate
    using Base.Threads
    using Parameters
    using SparseArrays
    using StaticArrays
    using Mmap


    # structures
    export FDMGrid, Signal1D
    # functions
    export gen3layv
    export rickerwave, signal1d, sourceposition, shotsposition
    export P2seis, seis2signals, imagecondition
    export discarray, todiscarray
    export propagate, propagateshots


    # standard constants
    const TAPER = 60
    const ATTENUATION = 0.0035  # found | Átila's = 0.008


    """
        FDMGrid(Δz::T, Δx::T, Δt::T, nz::Int, nx::Int, nt::Int)
    Structure for defining a 2D finite-difference grid.
    """
    struct FDMGrid{R, I}
        Δz::R
        Δx::R
        Δt::R
        nz::I
        nx::I
        nt::I
    end


    """
    Structure for 1D signal definition.
    """
    struct Signal1D{T}
        z::Integer
        x::Integer
        t1::Integer
        values::AbstractArray{T, 1}
    end


    """
        signal1d(z, x, values, t1=1)
    Function to conveniently create Signal1D structs.
    """
    signal1d(z, x, values, t1=1) = Signal1D(z, x, t1, values)


    """
        offset_signal(signal::Signal1D, offset::Integer)
    Offsets a Signal1D struct's position by some ammount.
    """
    function offset_signal(signal::Signal1D, offset::Integer)
        _signal = Signal1D(signal.z+offset,
                           signal.x+offset,
                           signal.t1,
                           signal.values)
    end

    """
        offset_signals(signals, offset::Integer)
    Offsets one or an array of Signal1D's some ammount. The return type is
    always Array{Signal1D, N}
    """
    function offset_signals(signals, offset)
        if typeof(signals) <: Signal1D
            [offset_signal(signals, offset)]
        else
            offset_signal.(signals, offset)
        end
    end


    """
        discarray(io, mode::String="r", type::DataType=Float64, dims=())
    This function is an interface to memory map creation. This interface has a
    header of Int64 values, the first one of them being the number of dimensions
    of the described array, _ndims, and then its dimensions. The body comes in
    sequence, and is a flattened version of that array.
    """
    function discarray(io, mode::String="r", type::DataType=Float64, dims=();
                       shared::Bool=false, pos::Integer=position(io))
        pos !== position(io) && seek(io, pos)
        if occursin("w", mode)
            _dims = Tuple(Int64(d) for d in dims)
            _ndims = Int64(length(_dims))
            write(io, _ndims, _dims...)
        else
            _ndims = read(io, Int64)
            _dims = Tuple(read(io, Int64) for i in 1:_ndims)
        end
        A = Mmap.mmap(io, Array{type, _ndims}, _dims)
        return(A)
    end

    function discarray(filename::String, mode::String="r",
                       type::DataType=Float64, dims=();
                       shared::Bool=false, pos::Integer=0, closeio::Bool=true)
        occursin("w", mode) && @assert dims !== ()
        io = open(filename, mode)
        A = discarray(io, mode, type, dims; shared=shared)
        closeio && close(io)
        return(A)
    end



    """
        todiscarray(filename::String, A::AbstractArray)
    This function takes an array, A, and creates a memory map with read and
    write permissions at a file. This file contains A elements as content after
    discarray header.
    """
    function todiscarray(io, A::AbstractArray)
        discarray(io, "w+", eltype(A), size(A)) .= A
    end


    """
        sourceposition(array::String, nz::Integer, nx::Integer)
    Get's source position as a tuple of elements (sz, sx), according to the
    seismic array gotten as parammeter. Possible array values are:
        - "endon",
        - "split",
        - "center".
    """
    function sourceposition(array::String, z, x)
        array==="endon"  && (pos=(z[1],          x[1]))
        array==="split"  && (pos=(z[1],          reduce(+,x)÷2))
        array==="center" && (pos=(reduce(+,x)÷2, reduce(+,x)÷2))
        return pos
    end

    sourceposition(array::String, x) = sourceposition(array, 1, x)

    function receptorpositions(array::String, n, Δx, sx)
        if array==="endon"
            span = collect(sx+Δx : Δx : sx+n*Δx)
        elseif array==="split"
            span = [collect(sx-(n÷2)*Δx : Δx : sx-Δx)
                    collect(sx+Δx : Δx : sx+(n-n÷2)*Δx)]
        elseif array==="center"
            receptorpositions(array, n, Δx, sx)
        end
    end



    """
        norm_ricker(t, σ)
    Computes the negative normalized ricker function at t with standard
    deviation σ.
    """
    ricker(t, σ) = 2/(√(3σ)π^(1/4)) * (1-(t/σ)^2) * exp(-t^2/(2σ^2))


    """
        ricker(t)
    Computes the negative unnormalized ricker function at t.
    """
    ricker(t) = (1-2t^2) * exp(-t^2)


    """
        rickerwave(ν::Real, Δt::Real)
    Compute a ricker wave with main frequency ν, sampled in time intervals Δt.
    For a stable signal assert ``ν << \\frac{1}{2 Δt}``.
    """
    function rickerwave(ν::Real, Δt::Real)
        @assert ν < 0.2 * 1.0/(2.0*Δt)
        interval = 2 * (2.2/(ν*Δt)) ÷ 2
        t = -interval÷2:interval÷2 .* (π*ν*Δt)
        ricker.(t)
    end


    """
        padextremes(A::AbstractArray, padding::Integer)
    Extends extremes of A by a determined padding length.
    """
    function padextremes(A::AbstractArray, padding::Integer)
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
    padvalue(A::AbstractArray, padding::Integer=1, value=zero(eltype(A)))
    Pad A copy 's extremes with value.
    """
    function padvalue(A::AbstractArray, padding::Integer=0, value=zero(eltype(A)))

        _A = value .* ones(eltype(A), Tuple(dim+2padding for dim in size(A)))
        _A[1+padding:end-padding, 1+padding:end-padding, 1] .= A
        _A
    end

    """
        pad_zeros_add_axes(A::AbstractArray,
                           padding::Integer=1,
                           dim::Integer=1)
    Pad extremes of A and add extra axes to it.
    """
    function pad_zeros_add_axes(A::AbstractArray,
                                padding::Integer=1,
                                dim::Integer=1)
        _A = zeros(eltype(A), (size(A, 1)+2*padding, size(A, 2)+2*padding, dim))
        _A[1+padding:end-padding, 1+padding:end-padding, 1] .= A
        _A
    end


    function gen3layv(nz, nx, h1, h2, V1, V2, V3)
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


    function get_∇²_stencil(order)
        coefs = central_difference_coefs(2, order)
        stencil = SVector{length(coefs)}(coefs)
    end


    """
        ∇²(A, IL)
    Function get the laplacian of a point inside an array. In order to optimize code,
    it actually calculates the laplacian without actually dividing the values by h².
    Then it's a h² laplacian function.
    """
    function ∇²(A, I, stencil, Δz, Δx)
        z, x = Tuple(I)
        r = length(stencil) ÷ 2
        sumz = sumx = zero(eltype(A))
        @inbounds @simd for i in -r:r
            sumz += A[z+i, x] * stencil[r+i+1]
        end
        @inbounds @simd for i in -r:r
            sumx += A[z, x+i] * stencil[r+i+1]
        end
        sumz/Δz^2 + sumx/Δx^2
    end


    """
        new_p(IP, Iv, cur_P, old_P, v, Δz, Δx, Δt)
    Function that calculates the next P value at point IP in pressure field,
    considering propagation velocity at that equal to point v[Iv]; cur_P and
    old_P are, respectively, the current and old pressure field arrays.
    """
    function new_p(IP, Iv, cur_P, old_P, v, ∇²_stencil, Δz, Δx, Δt)
        2cur_P[IP] - old_P[IP] + v[Iv]^2 * Δt * ∇²(cur_P, IP, ∇²_stencil, Δz, Δx)
    end

    """
        new_p(IP, Iv, cur_P, old_P, v, Δz, Δx, Δt, attenuation_factor)
    Function that calculates the next P value at point IP in pressure field,
    considering propagation velocity at that point equal to v[Iv]. cur_P and
    old_P are, respectively, the current and old pressure field arrays. The
    result is multiplied by an attenuation_factor ∈ [0,1]. old_P is attenuated
    twice, as in the usual taper attenuation scheme.
    """
    function new_p(IP, Iv, cur_P, old_P, v, ∇²_stencil, Δz, Δx, Δt, attenuation_factor)
        attenuation_factor * begin
            2cur_P[IP]  +  v[Iv]^2 * Δt * ∇²(cur_P, IP, ∇²_stencil, Δz, Δx) -
            attenuation_factor * old_P[IP]
        end
    end


    function update_P!(new_P!, cur_P, old_P, _v,
                       ∇²_stencil, I∇²r, Δz, Δx, Δt,
                       attenuation_factors)
        @threads for Iv in CartesianIndices(_v)
            IP = Iv + I∇²r
            new_P![IP] = new_p(IP, Iv, cur_P, old_P, _v,
                               ∇²_stencil, Δz, Δx, Δt,
                               attenuation_factors[IP])
        end
    end


    """
        get_attenuation_factor(dist)
    Add description...
    """
    function get_attenuation_factor(depth, attenuation)
        exp(-(attenuation*depth)^2)
    end


    function get_attenuation_factors(A, taper, attenuation, offset=0)
        factors = ones(Float64, (size(A,1), size(A,2)))
        for i in 1:taper+offset
            depth = taper-i+1
            factors[i, i:end+1-i] .=
            factors[i:end+1-i, i] .=
            factors[end-i+1, i:end-i+1] .=
            factors[i:end-i+1, end-i+1] .=
                get_attenuation_factor(depth, attenuation)
        end
        factors
        # SMatrix{size(factors)...}(factors)  # doesn't work with parallelization?
    end


    function putsignals!(signals, t, cur_P!)
        @inbounds @simd for signal in signals
            t_signal = t - signal.t1  # for T in 2:nt
            if (signal.t1 <= t) && (t_signal <= length(signal.values))
                cur_P![signal.z, signal.x] = signal.values[t_signal]
            end
        end
    end


    function propagate(grid, v, signals,
                       P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                       Pfile="",
                       seisfile="",
                       stencilorder::Integer=2,
                       direct_only::Bool=false,)
        # TODO: add keyword to pre-load Pfile as mmap
        ntcache = 3
        onlyseis = seisfile !== ""

        @unpack Δz, Δx, Δt, nz, nx, nt = grid

        ∇²_stencil = get_∇²_stencil(stencilorder)
        ∇²r = length(∇²_stencil) ÷ 2

        padding = taper + ∇²r
        (_nz, _nx) = (n+2padding for n in (nz, nx))

        # defining pressure field
        _P = zeros(eltype(P0), (_nz, _nx, ntcache))
        # pressure field of interest
        P = @view _P[1+padding:end-padding, 1+padding:end-padding, :]
        P[:,:,1] .= P0

        # offseting signals
        _signals = offset_signals(signals, padding)

        # defining velocity field
        _v = padextremes(v, taper)
        if direct_only
            _v[:,taper+1:end-taper] .= repeat(v[1:1,:], nz+2taper, 1)
        end

        if onlyseis
            savedseis = discarray(seisfile, "w+", Float64, (nt, nx))
            savedseis[1,:] .= P[1,:,1]
        else
            savedP = discarray(Pfile, "w+", Float64, (nz, nx, nt))
        end

        attenuation_factors = get_attenuation_factors(_P, taper,
                                                      attenuation,
                                                      ∇²r)

        I∇²r = CartesianIndex(∇²r, ∇²r)
        @inbounds for T in eltype(nt)(2):nt
            old_t = mod1(T-2, ntcache)
            cur_t = mod1(T-1, ntcache)
            new_t = mod1(T,   ntcache)
            old_P = @view _P[:,:,old_t]
            cur_P = @view _P[:,:,cur_t]
            new_P = @view _P[:,:,new_t]

            # putting 1D signals
            putsignals!(_signals, T, cur_P)

            # solve P equation
            update_P!(new_P, cur_P, old_P, _v,
                      ∇²_stencil, I∇²r, Δz, Δx, Δt,
                      attenuation_factors)

            if (new_t === ntcache) | (T === nt)
                ntslice = (1:new_t) .+ ntcache*floor(Int64, (T-1)/ntcache)
                if onlyseis
                    savedseis[ntslice,:] .= P[1,:,1:new_t]'
                else
                    savedP[:,:,ntslice] .= P[:,:,1:new_t]
                end
            end
        end

        if onlyseis
            return savedseis
        else
            return savedP
        end
    end

    using PyPlot

    function propagateshots(grid, v, shotssignals, nrec=10, Δxrec=1,
                             P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                             Pfile,
                             seisfile,
                             multiseisfile,
                             array::String="split",
                             stencilorder::Integer=2)
        @unpack nt = grid
        nshots = length(shotssignals)
        nwavelets = nshots * nrec

        multiseis = discarray(multiseisfile, "w+", Float64, (nt, nwavelets))

        for (S, signals) in enumerate(shotssignals)
            recpositions = receptorpositions(array, nrec, Δxrec, signals[1].x)
            savedslice = nrec*(S-1) .+ (1:nrec)
            savedseis = @view multiseis[:,savedslice]
            seis = propagate(grid, v, signals, P0, taper, attenuation;
                             seisfile=seisfile, stencilorder=stencilorder, direct_only=false)
            print(size(seis))
            @async begin
                plt.imshow(seis)
                plt.show()
            end
            savedseis[:,:] .= seis[:,recpositions]
            plt.imshow(savedseis)
            plt.show()
        end
    end


    #function migrate(grid, v, shot_signals, nrec=10, Δxrec=1,
    #                 P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
    #                 Pfile,
    #                 directseisfile,
    #                 revPfile,
    #                 migratedfile,
    #                 multiseisfile="",
    #                 array::String="split",
    #                 stencilorder::Integer=2,
    #                 direct_only::Bool=false,)
    #end



    """
        P2seis(P)
    Slices a seismogram out of a pressure field, P{Any, 3}(nz, nx, nt),
    and returns a copy of it.
    """
    function P2seis(P)
        seis = copy(P[1,:,:]')
    end

    """
        seis2signals(seis)
    Transforms a seismogram array into an array of Signal1D's with reversed
    wavelets. This array is useful in future in RTM steps.
    """
    function seis2signals(seis)
        [Signal1D(1, x, 1, seis[end:-1:1, x]) for x in axes(seis, 2)]
    end

    """
        imagecondition(Pfile, reversed_Pfile, migrated_file)
    Perform image condition between to filenames pointing to binary files
    formated as discarray does.
    """
    function imagecondition(Pfile, reversed_Pfile, migrated_file)
        P = discarray(Pfile)
        reversed_P = discarray(reversed_Pfile)
        (nz, nx, nt) = size(P)
        migrated = discarray(migrated_file, "w+", Float64, (nz, nx))
        migrated .= 0.
        for t in 1:nt
            @views migrated[:,:] .+= P[:,:,t] .* reversed_P[:,:,nt-t+1]
        end
    end
end
