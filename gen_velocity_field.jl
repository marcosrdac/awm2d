using Images
include("src/acoustics2d.jl")
include("src/discarrays.jl")
include("parameters.jl")
using .Acoustics2D
using .Discarrays


V = [1500. 2500. 3500.]' # m/s
# V = [0. 3.]' # m/s


mode = "3lay"
# mode = "pic"


if mode === "pic"
    # inputfilename = "models/marmousifull.png"
    inputfilename = "models/ref.png"


    luma(rgb) = 0.2126rgb[1] + 0.7152rgb[2] + 0.0722rgb[3]

    function toluma(img)
        L = Array{Float64}(undef, (size(img,1), size(img,2)))
        for I in CartesianIndices(L)
            L[I] = luma(img[I,:])
        end
        return L
    end
    function img2arr(filename, vmin=0, vmax=1)
        img = load(filename)
        arr = permutedims(rawview(channelview(img)), (2,3,1))
        arr = putbetween(toluma(arr), vmin, vmax)
    end
    Vmin, Vmax = minimum(V), maximum(V)
    # darker is faster
    # v = img2arr(inputfilename, Vmax, Vmin)
    # lighter is faster
    v = img2arr(inputfilename, Vmin, Vmax)
elseif mode === "3lay"
    # NX, NZ = 321, 321
    NX, NZ = 300, 300
    H1 = H2 = NZ÷3
    V1, V2, V3 = V[1:3]
    v = gen3layv(NZ, NX, H1, H2, V1, V2, V3)
end


using PyPlot
plt.imshow(v)
plt.colorbar()
plt.show()


todiscarray(vfile, v)
