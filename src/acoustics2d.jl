module Acoustics2D
    using Base.Threads
    using SpecialFunctions
    using Parameters
    using StaticArrays
    include("discarrays.jl")
    using .Discarrays

    # structures
    export FDMGrid, Signal1D
    # functions
    export gen3layv, img2arr
    export rickerwave, sourceposition, receptorpositions, shotsposition
    export P2seis, seis2signals, imagecondition
    export zerotoone, putbetween
    export propagate, propagate_rem, propagateshots, migrate
    export surfaceindices, edgeindices
    export seisgain


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

    function zerotoone(A, type=Float64)
        _A = similar(A, Float64)
        _A .= A .- minimum(A)
        _A ./= maximum(_A)
    end

    function putbetween(A, vmin, vmax)
        (vmax-vmin) .* zerotoone(A) .+ vmin
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
        ∇²(A, I, stencil, Δz, Δx)
    2D laplacian of a point inside an array.
    """
    function ∇²(A, I, stencil, Δz=1, Δx=1)
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
        compute_∇²(A, stencil, Δz, Δx)
    Function get the laplacian field over an array.
    """
    function compute_∇²(A, stencil, Δz=1, Δx=1; out=similar(A))
        r = length(stencil) ÷ 2
        @threads for I in CartesianIndices((1+r:size(A,1)-r, 1+r:size(A,2)-r))
            out[I] = ∇²(A, I, stencil, Δz, Δx)
        end
        return out
    end


    # """
        # new_p(I, curP, oldP, v, Δz, Δx, Δt)
    # Function that calculates the next P value at point IP in pressure field,
    # considering propagation velocity at that equal to point v[Iv]; curP and
    # oldP are, respectively, the current and old pressure field arrays.
    # """
    # function new_p(I, curP, oldP, v, ∇²_stencil, Δz, Δx, Δt)
        # 2curP[I] - oldP[I] + (v[I]*Δt)^2 * ∇²(curP, I, ∇²_stencil, Δz, Δx)
    # end

    """
        new_p(IP, Iv, curP, oldP, v, Δz, Δx, Δt, attenuation_factor)
    Function that calculates the next P value at point IP in pressure field,
    considering propagation velocity at that point equal to v[Iv]. curP and
    oldP are, respectively, the current and old pressure field arrays. The
    result is multiplied by an attenuation_factor ∈ [0,1]. oldP is attenuated
    twice, as in the usual taper attenuation scheme.
    """
    function new_p(I, curP, oldP, v, ∇²_stencil, Δz, Δx, Δt, attenuation_factor)
        attenuation_factor * begin
            2curP[I]  +  v[I]^2 * Δt^2 * ∇²(curP, I, ∇²_stencil, Δz, Δx) - 
            oldP[I] * attenuation_factor
        end
    end


    function update_P!(newP!, curP, oldP, v,
                       ∇²_stencil, Δz, Δx, Δt,
                       attenuation_factors)
        r = length(∇²_stencil) ÷ 2
        @threads for I in CartesianIndices((1+r:size(v,1)-r, 1+r:size(v,2)-r))
            newP![I] = new_p(I, curP, oldP, v,
                             ∇²_stencil, Δz, Δx, Δt,
                             attenuation_factors[I])
        end
    end


    function Q(k, x)
        if     k>=4  2*(1+2x^2) * Q(k-2, x) - Q(k-4, x)
        elseif k==2  1+2x^2
        else         one(x)
        end
    end

    using PyPlot

    # when it works, add attenuation_factor...
    function update_P_rem!(newP!, P, PP, vel, padding,
                           coss, Q, QQ, lap,
                           ∇²_stencil, Δz, Δx,
                           R, M, J, attenuation_factors)
        coss .= 0
        for k = 1:M
            if k == 1
                @. QQ = P  #* attenuation_factors
            else
                compute_∇²(Q, ∇²_stencil, Δz, Δx; out=lap)
                if k == 2
                    @. QQ = Q + 2(vel/R)^2 * lap
                else
                    @. QQ = 2Q + 4(vel/R)^2 * lap - QQ
                end
            end
            @. coss += J[k] * QQ
            Q, QQ = QQ, Q
        end
        @. newP! = 2coss - PP  #* attenuation_factors
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
            @. factors[i, i:end+1-i] = factors[end-i+1, i:end-i+1] =
                get_attenuation_factor(depth, attenuation)
            @. factors[i:end+1-i, i] = factors[i:end-i+1, end-i+1] =
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
                # curP![signal.z, signal.x] = signal.values[signalT]
                curP![signal.z, signal.x] += signal.values[signalT]
            end
        end
    end

    function surfaceindices(slice)
        [CartesianIndex(1, x) for x in slice]
    end

    function edgeindices(nz, nx)
        [[CartesianIndex(z, x) for z in [1,nz] for x in 1:nx]
         [CartesianIndex(z, x) for z in 1:nz for x in [1,nx]]]
    end

    function propagate(grid, v, signals,
                       P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                       Pfile="",
                       seisfile="",
                       stencilorder::Integer=8,
                       recpositions=surfaceindices(1:grid.nx),
                       direct_only::Bool=false,)
        ntcache = 3
        onlyseis = seisfile !== ""

        @unpack Δz, Δx, Δt, nz, nx, nt = grid

        ∇²_stencil = get_∇²_stencil(stencilorder)
        ∇²r = length(∇²_stencil) ÷ 2

        padding = taper + ∇²r
        (_nz, _nx) = (n+2padding for n in (nz, nx))
        @show _nz, _nx
        @show nz, nx

        # defining pressure field
        _P = zeros(eltype(P0), (_nz, _nx, ntcache))
        # pressure field of interest
        P = @view _P[1+padding:end-padding, 1+padding:end-padding, :]
        P[:,:,1] .= P0

        # offseting signals
        _signals = offset_signals(signals, padding)

        # defining velocity field
        _v = padextremes(v, padding)
        if direct_only
            _v[:,padding+1:end-padding] .= repeat(v[1:1,:], nz+padding, 1)
        end

        if onlyseis
            nrec = length(recpositions)
            savedseis = discarray(seisfile, "w+", Float64, (nt, nrec))
            savedseis[1,:] .= P[recpositions, 1]
        else
            savedP = discarray(Pfile, "w+", Float64, (nz, nx, nt))
        end

        attenuation_factors = get_attenuation_factors(_P, padding,
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
                      ∇²_stencil, Δz, Δx, Δt,
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
    end


    function propagate_rem(grid, v, signals,
                           P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                           Pfile="",
                           seisfile="",
                           stencilorder::Integer=8,
                           recpositions=surfaceindices(1:grid.nx),
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
        _v = padextremes(v, padding)
        if direct_only
            _v[:,padding+1:end-padding] .= repeat(v[1:1,:], nz+padding, 1)
        end

        if onlyseis
            nrec = length(recpositions)
            savedseis = discarray(seisfile, "w+", Float64, (nt, nrec))
            savedseis[1,:] .= P[recpositions, 1]
        else
            savedP = discarray(Pfile, "w+", Float64, (nz, nx, nt))
        end

        attenuation_factors = get_attenuation_factors(_P, padding,
                                                      attenuation,
                                                      ∇²r)

        # REM part
        Vmax = maximum(v)
        R = π*Vmax*√(1/Δx^2 + 1/Δz^2)
        M = ceil(Int, R*Δt+1)  # M > R Δt
        J = [besselj(2k, R*Δt) for k in 0:M]

        @show Vmax
        @show R
        @show M
        @show J


        cosLΔt_P = zero(_P[:,:,1])
        curQ = similar(_P[:,:,1])
        newQ = similar(_P[:,:,1])
        ∇²curQ = similar(_P[:,:,1])

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
            update_P_rem!(newP, curP, oldP, _v, padding,
                          cosLΔt_P, curQ, newQ, ∇²curQ,
                          ∇²_stencil, Δz, Δx,
                          R, M, J,  # differ!
                          attenuation_factors)

            # update_P!(newP, curP, oldP, _v,
                      # ∇²_stencil, I∇²r, Δz, Δx, Δt,
                      # attenuation_factors)

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

    # using PyPlot
    function revpropagateimcond(grid, v, signals,
                                P0=zero(v), taper=TAPER, attenuation=ATTENUATION;
                                Pfile="",
                                image=zero(v),
                                stencilorder::Integer=8,
                                recpositions=surfaceindices(1:grid.nx),
                                direct_only::Bool=false,)

        @unpack Δz, Δx, Δt, nz, nx = grid
        P = discarray(Pfile, "r+", Float64)  # maybe r+ is needed
        nt = size(P, 3)

        ∇²_stencil = get_∇²_stencil(stencilorder)
        ∇²r = length(∇²_stencil) ÷ 2

        padding = taper + ∇²r
        (_nz, _nx) = (n+2padding for n in (nz, nx))

        # defining pressure field
        _revP = zeros(eltype(P0), (_nz, _nx, 3))
        # pressure field of interest
        revP = @view _revP[1+padding:end-padding, 1+padding:end-padding, :]
        revP[:,:,1] .= P0

        # offseting signals
        _signals = offset_signals(signals, padding)

        # defining velocity field
        _v = padextremes(v, padding)

        attenuation_factors = get_attenuation_factors(_revP, padding,
                                                      attenuation,
                                                      ∇²r)

        I∇²r = CartesianIndex(∇²r, ∇²r)
        @inbounds for T in Int64.(2:nt)
            oldt = mod1(T-2, 3)
            curt = mod1(T-1, 3)
            newt = mod1(T,   3)
            oldP = @view _revP[:,:,oldt]
            curP = @view _revP[:,:,curt]
            newP = @view _revP[:,:,newt]

            # putting 1D signals
            putsignals!(_signals, T, curP)

            # solve P equation
            update_P!(newP, curP, oldP, _v,
                      ∇²_stencil, Δz, Δx, Δt,
                      attenuation_factors)

            @threads for I in CartesianIndices(v)
                image[I] += revP[I,newt] * P[I,2+nt-T]  # init=2->2
            end

            # if T==300
                # plt.imshow(revP[:,:,newt]; aspect="auto")
                # plt.show()

                # plt.imshow(P[:,:,2+nt-T]; aspect="auto")
                # plt.show()

                # plt.imshow(image; aspect="auto")
                # plt.show()
            # end

        end
        return image
    end



    function propagateshots(grid, v, P0=zero(v);
                            seisfile,
                            multiseisfile,
                            shotssignals,
                            shotsrecpositions=[surfaceindices(1:grid.nx)
                                               for shot in shotssignals],
                            stencilorder::Integer=8,
                            taper=TAPER, attenuation=ATTENUATION)
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

    function migrate(grid, v, P0=zero(v);
                     multiseisfile,   # original data
                     shotssignals,    # each shot data
                     migratedfile,    # output migration product file
                     Pfile,           # propagation file
                     reversedPfile,   # seismograms retro propagated
                     shotsrecpositions=[surfaceindices(1:grid.nx)
                                        for shot in shotssignals],
                     stencilorder::Integer=8,
                     taper=TAPER, attenuation=ATTENUATION)

        multiseis = discarray(multiseisfile)
        migrated = discarray(migratedfile, "w+", Float64, size(v))

        pos = 1
        for (S, signals) in enumerate(shotssignals)

                # debug
                S |> println

            # propagating signals
            @time P = propagate(grid, v, signals, P0; Pfile=Pfile)

            # reverse-propagating seimograms
            recpositions = shotsrecpositions[S]
            nrecs = length(recpositions)
            seisslice = pos:pos+nrecs-1
            seis = @view multiseis[:,seisslice]
            seissignals = seis2signals(seis, recpositions)
            # revP = propagate(grid, v, seissignals; Pfile=reversedPfile)
            @time revpropagateimcond(grid, v, seissignals; Pfile=Pfile, image=migrated)

            pos += nrecs
        end
    end



    """
        P2seis(P)
    Slices a seismogram out of a pressure field, P{Any, 3}(nz, nx, nt),
    and returns a copy of it.
    """
    function P2seis(P, recpositions=surfaceindices(axes(P,2)))
        seis = copy(P[recpositions,:]')
    end

    """
        seis2signals(seis)
    Transforms a seismogram array into an array of Signal1D's with reversed
    wavelets. This array is useful in future in RTM steps.
    """
    function seis2signals(seis, recpositions=surfaceindices(axes(seis,2)))
        [Signal1D(pos[1], pos[2], seis[end:-1:1, x])
         for (x, pos) in enumerate(recpositions)]
    end

    """
    Perform image condition between to filenames pointing to binary files
    formated as discarray does.
    """
    function imagecondition!(P, revP, image!)
        _P = @view revP[:,:,end:-1:1]
        for T in axes(P,3)
            for I in CartesianIndices(axes(P)[1:2])
                Threads.@spawn image![I] += P[I,T] * _P[I,T]
            end
        end
        return image!
    end

    function imagecondition(P, revP)
        image=zeros(eltype(P), size(P)[1:2])
        imagecondition!(P, revP, image)
    end

    function cachedimagecondition!(P, revP, image!; ntcache)
        nt = size(P, 3)
        Tslices    = [T+ntcache<nt ? (T:T+ntcache-1) : (T:nt) 
                      for T in 1:ntcache:nt]
        revTslices = [T-ntcache>1 ? (T-ntcache:T-1) : (1:T)
                      for T in nt:-ntcache:1]
        for S in 1:length(Tslices)
            (Tslice, revTslice) = (Tslices[S], revTslices[S])

                # debug
                "copying arrays to memory..." |> println

            # try @view
            Pcopy    = deepcopy(P[:,:,Tslice])
            revPcopy = deepcopy(revP[:,:,revTslice])
            
                # debug
                @show S/length(Tslices)

            imagecondition!(Pcopy, revPcopy, image!)
        end
    end

end
