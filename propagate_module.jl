module Propagate
    using Base.Threads
    using Parameters
    using SparseArrays
    using StaticArrays
    using Mmap: mmap


    # structures
    export FDM_Grid, Signal1D
    # functions
    export gen_3lay_v
    export rickerwave, signal1d, sourceposition, shotsposition
    export P2seis, seis2signals, image_condition
    export discarray, todiscarray
    export propagate, propagate_save, propagate_save_seis, propagate_shots


    # standard constants
    const TAPER = 60
    const ATTENUATION = 0.0035  # found | Átila's = 0.008


    """
        FDM_Grid(Δz::T, Δx::T, Δt::T, nz::Int, nx::Int, nt::Int)
    Structure for defining a 2D finite-difference grid.
    """
    struct FDM_Grid{R, I}
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
        discarray(filename::String, mode::String="r",
                  type::DataType=Float64, dims=())
    This function is an interface to memory map creation. This interface has a
    header of Int64 values, the first one of them being the number of dimensions
    of the described array, _ndims, and then its dimensions. The body comes in
    sequence, and is a flattened version of that array.
    """
    function discarray(filename::String, mode::String="r", type::DataType=Float64, dims=())
        if occursin("w", mode)
            @assert dims !== ()
            @assert eltype(dims) <: Int64
            io = open(filename, mode)
            _ndims = length(dims)
            write(io, _ndims, dims...)
        else
            io = open(filename, mode)
            _ndims = read(io, Int64)
            dims = Tuple(read(io, Int64) for i in 1:_ndims)
        end
        A = mmap(io, Array{type, _ndims}, dims)
        close(io)
        return(A)
    end


    """
        todiscarray(filename::String, A::AbstractArray)
    This function takes an array, A, and creates a memory map with read and
    write permissions at a file. This file contains A elements as content after
    discarray header.
    """
    function todiscarray(filename::String, A::AbstractArray)
        _A = discarray(filename, "w+", eltype(A), size(A))
        _A .= A
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
            span = [collect(sx-Δx : -Δx : sx-(n÷2)*Δx)
                    collect(sx+Δx : +Δx : sx+n*Δx)]
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
        pad_extremes(A::AbstractArray, padding::Integer)
    Extends extremes of A by a determined padding length.
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

    for func in (:propagate,    :propagate_save,    :propagate_save_seis,)
        quote
    #==========================================================================#
    function $func(grid, v, signals, P0=zero(v),
                   taper=TAPER, attenuation=ATTENUATION;
                   filename::String,
                   stencil_order::Integer=2,
                   direct_only::Bool=false,)
        @unpack Δz, Δx, Δt, nz, nx, nt = grid

        ∇²_stencil = get_∇²_stencil(stencil_order)
        ∇²r = length(∇²_stencil) ÷ 2
        offset = taper + ∇²r

        _v = pad_extremes(v, taper)
        _signals = offset_signals(signals, offset)
        _P = pad_zeros_add_axes(P0, offset, 3)

        # pressure field of interest
        P = @view _P[1+offset:end-offset, 1+offset:end-offset, :]

        direct_only && @views _v[:,taper+1:end-taper] .= repeat(v[1:1,:],
                                                                nz + 2taper, 1)

        attenuation_factors = get_attenuation_factors(_P, taper,
                                                      attenuation,
                                                      ∇²r)

        @initialize_saved_arrays $func  # declares saved_P or saved_seis

        I∇²r = CartesianIndex(∇²r, ∇²r)
        @inbounds for T in 2:nt
            new_t = mod1(T,   3)
            cur_t = mod1(T-1, 3)
            old_t = mod1(T-2, 3)

            old_P = @view _P[:,:,old_t]
            cur_P = @view _P[:,:,cur_t]
            new_P = @view _P[:,:,new_t]

            # putting 1D signals
            @inbounds @simd for signal in _signals
                T_signal = T - signal.t1  # for T in 2:nt
                if (signal.t1 <= T) && (T_signal <= length(signal.values))
                    cur_P[signal.z, signal.x] = signal.values[T_signal]
                end
            end

            # solve P equation
            update_P!(new_P, cur_P, old_P, _v,
                      ∇²_stencil, I∇²r, Δz, Δx, Δt,
                      attenuation_factors)

            @savesnapifsave $func
        end
        @returnpropperarrayifnotsave $func
    end
    #==========================================================================#
        end |> eval
    end


    function propagate_shots(grid, v, shots_signals, nrec=10, Δxrec=1, P0=zero(v),
                             taper=TAPER, attenuation=ATTENUATION;
                             multi_seis_file::String,
                             P_file::String,
                             array::String="split",
                             stencil_order::Integer=2,
                             direct_only::Bool=false,)
        @unpack Δz, Δx, Δt, nz, nx, nt = grid

        nshots = length(shots_signals)
        nwavelets = nshots * nrec

        ∇²_stencil = get_∇²_stencil(stencil_order)
        ∇²r = length(∇²_stencil) ÷ 2
        offset = taper + ∇²r

        _shot_signals = []
        for signals in shots_signals
            push!(_shot_signals, offset_signals(signals, offset))
        end

        _v = pad_extremes(v, taper)
        _rev_P = pad_zeros_add_axes(P0, offset, 3)
        _P = discarray(P_file, "w+", Float64, dims=(size(_rev_P, 1),
                                                        size(_rev_P, 2),
                                                        nt))

        # pressure field of interest
        P = @view _P[1+offset:end-offset, 1+offset:end-offset, :]

        direct_only && @views _v[:,taper+1:end-taper] .= repeat(v[1:1,:],
                                                                nz + 2taper, 1)

        attenuation_factors = get_attenuation_factors(_P, taper,
                                                      attenuation,
                                                      ∇²r)

        saved_seis = discarray(multi_seis_file, "w+", Float64, (nt, nwavelets))
        # saved_seis[1,:] .= P[1,:,1]
        

        I∇²r = CartesianIndex(∇²r, ∇²r)
        for S in 1:nshots
            saved_seis_cur = @view saved_seis[1+(S-1)*nrec:S*nrec]
            recpositions = receptorpositions(array, n, Δxrec, _signals[1].x)
            _signals = _shot_signals[S]
            _P .=0
            P[:,:,0] .= P0

            for T in 2:nt
                new_t = mod1(T,   3)
                cur_t = mod1(T-1, 3)
                old_t = mod1(T-2, 3)

                old_P = @view _P[:,:,old_t]
                cur_P = @view _P[:,:,cur_t]
                new_P = @view _P[:,:,new_t]

                # putting 1D signals
                @inbounds @simd for signal in _signals
                    T_signal = T - signal.t1  # for T loops from 2, not 1
                    if (signal.t1 <= T) && (T_signal <= length(signal.values))
                        cur_P[signal.z, signal.x] = signal.values[T_signal]
                    end
                end

                # solve P equation
                update_P!(new_P, cur_P, old_P, _v,
                          ∇²_stencil, I∇²r, Δz, Δx, Δt,
                          attenuation_factors)

                saved_seis[T,recpositions] .= P[1,recpositions,new_t]
            end
        end
    end


    # function migrate(grid, v, signals, P0=zero(v),
                     # taper=TAPER, attenuation=ATTENUATION;
                     # P_file::String,
                     # rev_P_file::String,
                     # migrated_file::String,
                     # stencil_order::Integer=2,
                     # direct_only::Bool=false,)
        # @unpack Δz, Δx, Δt, nz, nx, nt = grid

        # ∇²_stencil = get_∇²_stencil(stencil_order)
        # ∇²r = length(∇²_stencil) ÷ 2
        # offset = taper + ∇²r

        # _v = pad_extremes(v, taper)
        # _v_direct = deepcopy(_v)
        # _v_direct[:,taper+1:end-taper] .= repeat(v[1:1,:])
        # _signals = offset_signals(signals, offset)
        # _P = pad_zeros_add_axes(P0, offset, 3)
        # # pressure field of interest
        # P = @view _P[1+offset:end-offset, 1+offset:end-offset, :]

        # attenuation_factors = get_attenuation_factors(_P, taper,
                                                      # attenuation,
                                                      # ∇²r)
    # end


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
        image_condition(P_file, reversed_P_file, migrated_file)
    Perform image condition between to filenames pointing to binary files
    formated as discarray does.
    """
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
