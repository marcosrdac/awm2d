include("./parameters.jl")

@time propagate_save(grid, P0, v, signal; filename=P_file)

using BenchmarkTools
#@btime propagate_save($grid, $P0, $v, $signal; filename=$P_file,)
#@btime propagate_save_seis($grid, $P0, $v, $signal; filename=$P_file,)
