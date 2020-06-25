using PyPlot

io = open("P.bin", "r")
all_shape = read(io, sizeof(Int)*3)
all_shape = reinterpret(Int, all_shape)
all_shape = Tuple(all_shape)
shape = (all_shape[1], all_shape[2])
println(shape)

read(io, 2500*sizeof(Float64)*shape[1]*shape[2])
A = read(io, sizeof(Float64)*shape[1]*shape[2])
A = reinterpret(Float64, A)
A = reshape(A, (shape[1], shape[2]))
close(io)

plt.imshow(A)
plt.show()

