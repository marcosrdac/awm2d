from os.path import basename
from sys import argv
import matplotlib.pyplot as plt
from src.discarrays import discarray

# defining parameters
nt = 1  # default
if len(argv) > 1:
    filename = argv[1]
    if len(argv) > 2:
        nt = int(argv[2])
else:
    filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/v.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/seis.bin"
    filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 200-1
    filename = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"

arr = discarray(filename, order='F')

# plotting data
if arr.ndim == 1:
    plt.plot(arr)
elif arr.ndim == 2:
    #plt.imshow(arr, aspect='auto')
    #plt.imshow(arr, aspect='auto', vmin=-.006, vmax=.006)
    #plt.imshow(arr, aspect='auto', vmin=-.01, vmax=.01)
    # plt.imshow(arr, aspect='auto', vmin=-.2, vmax=.2)
    plt.imshow(arr, aspect='auto', interpolation='nearest')
    plt.colorbar()
elif arr.ndim == 3:
    plt.imshow(arr[:, :, nt])
    plt.colorbar()


plt.title(basename(filename))
# plt.savefig("migrated.jpg")
plt.show()
