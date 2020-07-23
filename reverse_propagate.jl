using InteractiveUtils
include("./reff.jl")
include("./parameters.jl")

# 2D signal parameters
begin
    P = open_mmap(P_file)
    seis = P[1,:,:]'
    position = CartesianIndex(1, 1)
    signal = Signal2D(seis, position)
end

@time propagate(grid, P0, v, signal;
                save=true,
                filename=reversed_P_file,
                only_seis=false)
