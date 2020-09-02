import numpy as np
from os.path import basename
import matplotlib.pyplot as plt

# defining parameters
filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
filename = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
filename = "/mnt/hdd/home/tmp/awp_data/v.bin"
filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 3000-1
filename = "/mnt/hdd/home/tmp/awp_data/migrated.bin"

# reading file array
io = open(filename, "r")
n = np.fromfile(io, dtype=np.int64, count=1)[0]
shape = np.fromfile(io, dtype=np.int64, count=n)

if len(shape) == 2:
    arr = np.fromfile(io,
                       dtype=np.float64,
                       count=shape[0]*shape[1]*np.float64().itemsize)
    arr = arr.reshape((shape[0], shape[1]), order='F')
elif len(shape) == 3:
    arr = np.fromfile(io,
                    dtype=np.float64,
                    count=shape[0]*shape[1],
                    offset=int(shape[0]*shape[1]*nt)*np.float64().itemsize)
    arr = arr.reshape((shape[0], shape[1]), order='F')
else:
    raise Exception("Shape not understood.")


# plotting data
if len(shape) == 2:
    #plt.imshow(arr, aspect='auto')
    #plt.imshow(arr, aspect='auto', vmin=-.006, vmax=.006)
    #plt.imshow(arr, aspect='auto', vmin=-.01, vmax=.01)
    #plt.imshow(arr, aspect='auto', vmin=-.2, vmax=.2)
    plt.imshow(arr, aspect='auto')
else:
    plt.imshow(arr)


plt.colorbar()
plt.title(basename(filename))
plt.savefig("migrated.jpg")
plt.show()
