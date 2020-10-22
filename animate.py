#!/usr/bin/env python3

from os.path import join, basename, splitext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from src.discarrays import discarray


def animate_snaps(in_fn, out_fn, order='C',
                  norm="frameabsmax",
                  timestep=30, interval=60,
                  cmap='seismic'):
    '''
    Saves a [nz,nx,nt] array time snaps as a gif.
    Norm can be:
        - "fullminmax";
        - "frameminmax";
        - "fullabsmax";
        - "frameabsmax";
        - (vmin, vmax);
    '''
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_axes([0, 0, 1, 1])

    arr = discarray(in_fn, "r", order=order)
    nz, nx, nt = arr.shape
    frame = arr[:, :, 0]

    im = ax.imshow(frame, aspect='auto', cmap=cmap)

    vmin = vmax = None
    if norm is not None:
        if "full" in norm:
            vmin, vmax = np.min(arr), np.max(arr)
        elif "frame" in norm:
            vmin, vmax = np.min(frame), np.max(frame)
        elif isinstance(norm, tuple):
            vmin, vmax = norm

        if "absmax" in norm:
            absmax = np.max(np.abs([vmin, vmax]))
            vmin, vmax = -absmax, absmax
        im.set_clim(vmin, vmax)

    def update(t, vmin, vmax):
        frame = arr[:, :, t]
        im.set_data(frame)
        if norm is not None:
            if "frame" in norm:
                vmin, vmax = np.min(frame), np.max(frame)
                if "absmax" in norm:
                    absmax = np.max(np.abs([vmin, vmax]))
                    vmin, vmax = -absmax, absmax
            im.set_clim(vmin=vmin, vmax=vmax)
        return(ax)

    anim = FuncAnimation(fig, update, frames=np.arange(1, nt, timestep),
                         interval=interval, fargs=(vmin, vmax))
    anim.save(out_fn, dpi=80, writer='imagemagick')
    # plt.show()


if __name__ == '__main__':
    folder = "/mnt/hdd/home/tmp/awp_data"
    Pfile = join(folder, "P.bin")
    revPfile = join(folder, "reversed_P.bin")

    outfolder = "animations"

    # files = [Pfile, revPfile]
    files = [Pfile]
    # files = [revPfile]

    for f in files:
        outf = join(outfolder, f"{splitext(basename(f))[0]}.gif")
        print(outf)
        animate_snaps(f, outf, order='F')
