module Discarrays
    using Mmap

    export discarray, todiscarray

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

end
