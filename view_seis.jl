using PyPlot

io = open("seis.bin", "r")
shape = read(io, 2*sizeof(Int))
shape = reinterpret(Int, shape)
shape = Tuple(shape)

A = read(io, sizeof(Float64)*shape[1]*shape[2])
A = reinterpret(Float64, A)
A = reshape(A, shape)
close(io)

#plt.imshow(A[1300:end, :]; aspect="auto", vmax=.1)
plt.imshow(A[:, :]; aspect="auto", vmax=.1)
plt.colorbar()
plt.show()
