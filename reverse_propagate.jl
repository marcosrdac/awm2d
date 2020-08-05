using Base.Threads
include("./parameters.jl")

# 2D signal parameters
begin
    P = open_mmap(P_file)
    seis = P[1,:,:]'
    direct_seis = open_mmap(direct_seis_file)
    seis_wo_direct = similar(seis)
    seis_wo_direct .= seis
    @threads for I in eachindex(direct_seis)
        seis_wo_direct[I] -= direct_seis[I]
    end
    
    position = CartesianIndex(1, 1)
    #signal = Signal2D(seis, position)
    signal = Signal2D(seis_wo_direct, position)
end

@time propagate_2d_signal(grid, P0, v, signal;
                          save=true,
                          filename=reversed_P_file,
                          only_seis=false)
