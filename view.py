#!/usr/bin/env python3


from os.path import basename
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from src.discarrays import discarray
from src.seisgain import seisgain

# defining parameters
nt = 1  # default
if len(argv) > 1:
    filename = argv[1]
    if len(argv) > 2:
        nt = int(argv[2])
else:
    filename = "/mnt/hdd/home/tmp/awp_data/model.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/direct_seis.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/source_signature.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/3lay_migrated.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/migrated_marmousi_5.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/v.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/migrated.bin"
    filename = "/mnt/hdd/home/tmp/awp_data/multi_seis.bin"
    # filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 2600-1
    filename, nt = "/mnt/hdd/home/tmp/awp_data/P.bin", 500-1
    filename = "/mnt/hdd/home/tmp/awp_data/seis.bin"

arr = discarray(filename, order='F')

# plotting data
if arr.ndim == 1:
    plt.plot(arr)
elif arr.ndim == 2:
    # plt.imshow(arr, aspect='auto')
    #plt.imshow(arr, aspect='auto', vmin=-.006, vmax=.006)
    #plt.imshow(arr, aspect='auto', vmin=-.01, vmax=.01)
    # plt.imshow(arr, aspect='auto', vmin=-.2, vmax=.2)
    # plt.imshow(arr, aspect='auto', interpolation='nearest')
    # plt.imshow(arr, aspect='auto', cmap='seismic', interpolation='nearest')
    # plt.imshow(arr, aspect='auto', cmap='seismic',)

    # plt.imshow(arr, aspect='auto', cmap='seismic', vmin=-0.06, vmax=0.06)
    # plt.imshow(arr, aspect='auto', cmap='viridis', vmin=-0.06, vmax=0.06)
    # plt.imshow(seisgain(arr, a=2.), aspect='auto', cmap='seismic')
    # plt.imshow(arr, aspect='auto')

    # _arr = seisgain(arr, a=1.8, b=-1.)
    _arr = seisgain(arr, a=2., b=0.)
    # _arr = arr

    # plt.subplot(211)
    # plt.imshow(discarray("/mnt/hdd/home/tmp/awp_data/v.bin", order='F'), cmap='Greys_r')
    # plt.colorbar()
    # plt.subplot(212)
    # plt.imshow(_arr, cmap='Greys_r')

    plt.imshow(_arr, aspect='auto', cmap='Greys_r')
    plt.colorbar()
elif arr.ndim == 3:
    plt.imshow(arr[:, :, nt])
    plt.colorbar()


plt.title(basename(filename))
plt.tight_layout()
plt.savefig("savedfigs/saved.jpg")
plt.show()
