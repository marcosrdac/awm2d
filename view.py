from os.path import basename
import matplotlib.pyplot as plt
from discarray import discarray


# defining parameters
filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
filename = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
filename = "/mnt/hdd/home/tmp/awp_data/v.bin"
filename = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
filename = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"
filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 1-1

arr = discarray(filename, order='F')

# plotting data
if arr.ndim == 1:
    plt.plot(arr)
elif arr.ndim == 2:
    #plt.imshow(arr, aspect='auto')
    #plt.imshow(arr, aspect='auto', vmin=-.006, vmax=.006)
    #plt.imshow(arr, aspect='auto', vmin=-.01, vmax=.01)
    #plt.imshow(arr, aspect='auto', vmin=-.2, vmax=.2)
    plt.imshow(arr, aspect='auto')
elif arr.ndim == 3:
    plt.imshow(arr[:, :, nt])


plt.colorbar()
plt.title(basename(filename))
# plt.savefig("migrated.jpg")
plt.show()
