include("src/acoustics2d.jl")
include("src/discarrays.jl")
include("parameters.jl")
using .Acoustics2D
using .Discarrays


# V = [3. 5. 9.]' # km/s
V = [0. 3.]' # km/s

mode = "pic"


if mode === "pic"
    # inputfilename = "models/marmousifull.png"
    inputfilename = "models/dif_model.png"
    Vmin, Vmax = minimum(V), maximum(V)
    # darker is faster
    # v = img2arr(inputfilename, Vmax, Vmin)
    # darker is slower
    v = img2arr(inputfilename, Vmin, Vmax)
elseif mode === "3lay"
    NX, NZ = 321, 321
    # NX, NZ = 50, 50
    H1 = H2 = NZ÷3
    V1, V2, V3 = V[1:3]
    v = gen3layv(NZ, NX, H1, H2, V1, V2, V3)
end


using PyPlot
plt.imshow(v)
plt.colorbar()
plt.show()


todiscarray(vfile, v)
