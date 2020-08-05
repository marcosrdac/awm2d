using Base.Threads
include("./parameters.jl")

begin
    P = open_mmap(P_file)
    seis = slice_seismogram(P)
    direct_seis = open_mmap(direct_seis_file)

    # removing direct_wave from seismogram
    @threads for I in eachindex(direct_seis)
        seis[I] -= direct_seis[I]
    end
    
    signal2d = Signal2D(seis, position, CartesianIndex(1, 1))
end

@time propagate_2d(grid, P0, v, signal;
                   filename=reversed_P_file,)
