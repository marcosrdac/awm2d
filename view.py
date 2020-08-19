import numpy as np
from os.path import basename
import matplotlib.pyplot as plt

# defining parameters
filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
filename = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
filename = "/mnt/hdd/home/tmp/awp_data/v.bin"
filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 10
filename = "/mnt/hdd/home/tmp/awp_data/migrated.bin"

# reading file array
io = open(filename, "r")
n = np.fromfile(io, dtype=np.int64, count=1)[0]
shape = np.fromfile(io, dtype=np.int64, count=n)

if len(shape) == 2:
    seis = np.fromfile(io,
                       dtype=np.float64,
                       count=shape[0]*shape[1]*np.float64().itemsize)
    seis = seis.reshape((shape[0], shape[1]), order='F')
else:
    P = np.fromfile(io,
                    dtype=np.float64,
                    count=shape[0]*shape[1],
                    offset=int(shape[0]*shape[1]*nt)*np.float64().itemsize)
    P = P.reshape((shape[0], shape[1]), order='F')


# plotting data
if len(shape) == 2:
    #plt.imshow(seis, aspect='auto')
    #plt.imshow(seis, aspect='auto', vmin=-.006, vmax=.006)
    #plt.imshow(seis, aspect='auto', vmin=-.01, vmax=.01)
    #plt.imshow(seis, aspect='auto', vmin=-.2, vmax=.2)
    plt.imshow(seis, aspect='auto')
else:
    plt.imshow(P)


plt.colorbar()
plt.title(basename(filename))
# plt.savefig("migrated.jpg")
plt.show()
