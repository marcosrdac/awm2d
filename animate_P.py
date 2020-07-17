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

    def update(t):
        cur_P = np.fromfile(io,
                            dtype=np.float64,
                            count=nz*nx,
                            offset=int(nz*nx*t) * np.float64().itemsize)
        cur_P = cur_P.reshape((nz, nx), order='F')
        im.set_data(cur_P)
        return(ax)

    anim = FuncAnimation(fig, update, frames=np.arange(0, nt, 10), interval=70)
    anim.save(out_fn, dpi=80, writer='imagemagick')
    plt.show()


if __name__ == '__main__':
    gen_P_gif("data/P.bin", "P.gif")
