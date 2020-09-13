module Acoustics2D
    using Base.Threads
    using Parameters
    using StaticArrays
    using Images
    include("discarrays.jl")
    using .Discarrays

    # structures
    export FDMGrid, Signal1D
    # functions
    export gen3layv, img2arr
    export rickerwave, sourceposition, receptorpositions, shotsposition
    export P2seis, seis2signals, imagecondition
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
        startT::Integer
        values::AbstractArray{T, 1}
    end

    Signal1D(z, x, values) = Signal1D(z, x, 1, values)


    """
        offset_signal(signal::Signal1D, offset::Integer)
    Offsets a Signal1D struct's position by some ammount.
    """
    function offset_signal(signal::Signal1D, offset::Integer)
        _signal = Signal1D(signal.z+offset,
                           signal.x+offset,
                           signal.startT,
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
        sourceposition(array::String, nz::Integer, nx::Integer)
    Get's source position as a tuple of elements (sz, sx), according to the
    seismic array gotten as parammeter. This function is useful for one-shot
    propagation products. Possible array values are:
        - "endon" or "endon2right";
        - "endon2left";
        - "split"; and
        - "center".
    """
    function sourceposition(array::String, z, x)
        array==="endon"       && (pos=(z[1],          x[1]))
        array==="endon2right" && (pos=(z[1],          x[1]))
        array==="endon2left"  && (pos=(z[1],          x[end]))
        array==="split"       && (pos=(z[1],          reduce(+,x)÷2))
        array==="center"      && (pos=(reduce(+,x)÷2, reduce(+,x)÷2))
        return pos
    end

    sourceposition(array::String, x) = sourceposition(array, 1, x)

    """
        receptorpositions(array::String, nz::Integer, nx::Integer)
    Get's valid receptor positions as an array of CartesianIndex elements
    (sz, sx), according to the seismic array gotten as parammeter. Possible 
    array values are:
        - "endon" or "endon2right";
        - "endon2left";
        - "split"; and
        - "center".
    TODO: receptors directly above source, if source is buried (if it makes
    any sense to implement that...)
    """
    function receptorpositions(array::String, sx, Δx, n, nx=Inf)
        if (array === "endon") | (array === "endon2right")
            span = [CartesianIndex(1, x)
                    for x in sx+Δx : Δx : sx+n*Δx
                    if 1 <= x <= nx]
        elseif array === "endon2left"
            span = [CartesianIndex(1, x)
                    for x in sx-n*Δx : Δx : sx-Δx
                    if 1 <= x <= nx]
        elseif array === "split"
            span = [
                    [CartesianIndex(1, x)
                     for x in sx-(n÷2)*Δx : Δx : sx-Δx
                     if 1 <= x <= nx]
                    [CartesianIndex(1, x)
                     for x in sx+Δx : Δx : sx+(n-n÷2)*Δx
                     if 1 <= x <= nx]
                   ]
        elseif array === "center"
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
        len = 2 * (2.2/(ν*Δt)) ÷ 2
        t = (π*ν*Δt) .* (-len÷2:len÷2)
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


    luma(rgb) = 0.2126rgb[1] + 0.7152rgb[2] + 0.0722rgb[3]

    function toluma(img)
        L = Array{Float64}(undef, (size(img,1), size(img,2)))
        for I in CartesianIndices(L)
            L[I] = luma(img[I,:])
        end
        return L
    end

    function zerotoone(A, type=Float64)
        _A = similar(A, Float64)
        _A .= A .- minimum(A)
        _A ./= maximum(_A)
    end

    function putbetween(A, vmin, vmax)
        (vmax-vmin) .* zerotoone(A) .+ vmin
    end

    function img2arr(filename, vmin=0, vmax=1)
        img = load(filename)
        arr = permutedims(rawview(channelview(img)), (2,3,1))
        arr = putbetween(toluma(arr), vmin, vmax)
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
        new_p(IP, Iv, curP, oldP, v, Δz, Δx, Δt)
    Function that calculates the next P value at point IP in pressure field,
    considering propagation velocity at that equal to point v[Iv]; curP and
    oldP are, respectively, the current and old pressure field arrays.
    """
    function new_p(IP, Iv, curP, oldP, v, ∇²_stencil, Δz, Δx, Δt)
        2curP[IP] - oldP[IP] + v[Iv]^2 * Δt * ∇²(curP, IP, ∇²_stencil, Δz, Δx)
    end

    """
        new_p(IP, Iv, curP, oldP, v, Δz, Δx, Δt, attenuation_factor)
    Function that calculates the next P value at point IP in pressure field,
    considering propagation velocity at that point equal to v[Iv]. curP and
    oldP are, respectively, the current and old pressure field arrays. The
    result is multiplied by an attenuation_factor ∈ [0,1]. oldP is attenuated
    twice, as in the usual taper attenuation scheme.
    """
    function new_p(IP, Iv, curP, oldP, v, ∇²_stencil, Δz, Δx, Δt, attenuation_factor)
        attenuation_factor * begin
            2curP[IP]  +  v[Iv]^2 * Δt * ∇²(curP, IP, ∇²_stencil, Δz, Δx) -
            attenuation_factor * oldP[IP]
        end
    end


    function update_P!(newP!, curP, oldP, _v,
                       ∇²_stencil, I∇²r, Δz, Δx, Δt,
                       attenuation_factors)
        @threads for Iv in CartesianIndices(_v)
            IP = Iv + I∇²r
            newP![IP] = new_p(IP, Iv, curP, oldP, _v,
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
            factors[i, i:end+1-i] .= factors[end-i+1, i:end-i+1] .=
                get_attenuation_factor(depth, attenuation)
            factors[i:end+1-i, i] .= factors[i:end-i+1, end-i+1] .=
                get_attenuation_factor(depth, attenuation)
        end
        factors
        # SMatrix{size(factors)...}(factors)  # doesn't work with parallelization?
    end


    function putsignals!(signals, newT, curP!, startT=2)
        @inbounds @simd for signal in signals
            curT = one(newT) + newT - startT
            signalT = one(newT) + curT - signal.startT
            if (signal.startT <= curT) & (signalT <= length(signal.values))
                curP![signal.z, signal.x] = signal.values[signalT]
            end
        end
    end
    
    function surfaceindexes(slice)
        [CartesianIndex(1, x) for x in slice]
    end

    function propagate(grid, v, signals,
                       P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                       Pfile="",
                       seisfile="",
                       stencilorder::Integer=8,
                       recpositions=surfaceindexes(1:grid.nx),
                       direct_only::Bool=false,)
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
            nrec = length(recpositions)
            savedseis = discarray(seisfile, "w+", Float64, (nt, nrec))
            savedseis[1,:] .= P[recpositions, 1]
        else
            savedP = discarray(Pfile, "w+", Float64, (nz, nx, nt))
        end

        attenuation_factors = get_attenuation_factors(_P, taper,
                                                      attenuation,
                                                      ∇²r)

        I∇²r = CartesianIndex(∇²r, ∇²r)
        @inbounds for T in eltype(nt)(2):nt
            oldt = mod1(T-2, ntcache)
            curt = mod1(T-1, ntcache)
            newt = mod1(T,   ntcache)
            oldP = @view _P[:,:,oldt]
            curP = @view _P[:,:,curt]
            newP = @view _P[:,:,newt]

            # putting 1D signals
            putsignals!(_signals, T, curP)

            # solve P equation
            update_P!(newP, curP, oldP, _v,
                      ∇²_stencil, I∇²r, Δz, Δx, Δt,
                      attenuation_factors)

            if (newt === ntcache) | (T === nt)
                ntslice = (1:newt) .+ ntcache*floor(Int64, (T-1)/ntcache)
                if onlyseis
                    savedseis[ntslice,:] .= P[recpositions,1:newt]'
                else
                    savedP[:,:,ntslice] .= P[:,:,1:newt]
                end
            end
        end

        if onlyseis
            return savedseis
        else
            return savedP
        end
    end

    function propagateshots(grid, v, shotssignals,
                            shotsrecpositions=[surfaceindexes(1:grid.nx)
                                               for shot in shotssignals],
                            P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                            seisfile,
                            multiseisfile,
                            stencilorder::Integer=8)
        @unpack nt = grid

        # nshots = length(shotssignals)
        nwavelets = reduce(+, [length(recs) for recs in shotsrecpositions])

        multiseis = discarray(multiseisfile, "w+", Float64, (nt, nwavelets))

        pos = 1
        for (S, signals) in enumerate(shotssignals)
            signals = shotssignals[S]
            recpositions = shotsrecpositions[S]
            
            seis = propagate(grid, v, signals, P0, taper, attenuation;
                             seisfile=seisfile, stencilorder=stencilorder,
                             recpositions=recpositions)

            nrecs = length(recpositions)
            savedslice = pos:pos+nrecs-1

            multiseis[:,savedslice] .= seis

            pos += nrecs
        end
    end


    function migrate(grid, v, shotssignals,
                     shotsrecpositions=[surfaceindexes(1:grid.nx)
                                        for shot in shotssignals],
                     P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                     Pfile,           # signal propagated
                     directseisfile,  # direct signal propagated
                     multiseisfile,   # seismograms
                     revPfile,        # seismograms retro propagated
                     migratedfile,    # final migration product file
                     stencilorder::Integer=8)
    end



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
        [Signal1D(1, x, seis[end:-1:1, x]) for x in axes(seis, 2)]
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
