import numpy as np
import matplotlib.pyplot as plt

f = open("P.bin", "r")
shape = np.fromfile(f, dtype=np.uint64, count=3)
for i in range(shape[2]):
    arr = np.fromfile(f, dtype=np.float64, count=shape[0]*shape[1])
    arr = arr.reshape((shape[1], shape[0])).T
    plt.imshow(arr)
    plt.colorbar()
    plt.show()
