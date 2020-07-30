#!/usr/bin/env python3

from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from matplotlib.animation import FuncAnimation


def gen_P_gif(in_fn, out_fn):
    '''
    Saves a [nz,nx,nt] array time snaps as a gif.
    '''
    # opening read-only file as IO
    io = open(in_fn, "r")
    # read number of dimensions
    n = np.fromfile(io, dtype=np.int64, count=1)[0]
    # read dimensions as shape
    shape = np.fromfile(io, dtype=np.int64, count=n)
    nz, nx, nt = shape

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])

    cur_P = np.fromfile(io, dtype=np.float64, count=nz*nx)
    cur_P = cur_P.reshape((nz, nx), order='F')

    im = ax.imshow(cur_P, aspect='auto', vmin=-.1, vmax=.1)

    time_sep = 10

    def update(t):
        cur_P = np.fromfile(io,
                            dtype=np.float64,
                            count=nz*nx,
                            offset=int(nz*nx*time_sep) * np.float64().itemsize)
        try:
            cur_P = cur_P.reshape((nz, nx), order='F')
            im.set_data(cur_P)
        except:
            print(nt)
            return()
        return(ax)

    anim = FuncAnimation(fig, update, frames=np.arange(1, nt-1, time_sep), interval=70)
    anim.save(out_fn, dpi=80, writer='imagemagick')
    #plt.show()


if __name__ == '__main__':
    P_file = "/mnt/hdd/home/tmp/awp_data/P.bin"
    P_sub_direct_file = "/mnt/hdd/home/tmp/awp_data/P_sub_direct.bin"
    reversed_P_file = "/mnt/hdd/home/tmp/awp_data/reversed_P.bin"

    #gen_P_gif(P_file,            "animations/P.gif")
    gen_P_gif(P_sub_direct_file, "animations/P_sub_direct.gif")
    #gen_P_gif(reversed_P_file,   "animations/reversed_P.gif")
