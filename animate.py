#!/usr/bin/env python3

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from src.discarrays import discarray


def animate_snaps(in_fn, out_fn, order='C'):
    '''
    Saves a [ny,nx,nt] array time snaps as a gif.
    '''
    # opening read-only file as IO
    arr = discarray(in_fn, "r", order=order)
    nz, nx, nt = arr.shape

    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_axes([0, 0, 1, 1])

    cur_P = arr[:, :, 0]
    im = ax.imshow(cur_P, aspect='auto', vmin=-.1, vmax=.1)

    time_sep = 30

    def update(t):
        cur_P = arr[:, :, t]
        im.set_data(cur_P)
        # NORMALIZATION METHODS
        #   1. NONE: comment every other method
        #
        # min, max: useful for methods 2 and 3
        vmin, vmax = np.min(cur_P), np.max(cur_P)
        #   2. min=min, max=max
        # im.set_clim(vmin=vmin, vmax=vmax)
        #
        #   3. min=-maxabs, max=+maxabs
        maxabs = np.max(np.abs([vmin, vmax]))
        im.set_clim(vmin=-maxabs, vmax=maxabs)
        return(ax)

    anim = FuncAnimation(fig, update, frames=np.arange(
        1, nt, time_sep), interval=60)
    anim.save(out_fn, dpi=80, writer='imagemagick')
    # plt.show()


if __name__ == '__main__':
    folder = "/mnt/hdd/home/tmp/awp_data"

    animate_snaps(join(folder, "P.bin"),
                  "animations/P.gif", order='F')

    # animate_snaps(join(folder, "reversed_P.bin"),
    #              "animations/reversed_P.gif", order='F')
