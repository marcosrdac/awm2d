using InteractiveUtils
include("./reff.jl")
include("./parameters.jl")

@time image_condition(P_file, reversed_P_file, migrated_file)
# 102 s
# 149 s
