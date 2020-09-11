include("./parameters.jl")
include("./propagate_module.jl")

using .Propagate

# three layered model parameters
H1 = (1*NZ)÷3
H2 = (1*NZ)÷3
V1 = 3. # km/s
V2 = 5. # km/s
V3 = 9. # km/s

v = gen3layv(NZ, NX, H1, H2, V1, V2, V3)
vdisc = todiscarray(vfile, v)
