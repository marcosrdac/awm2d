#!/usr/bin/env python3

from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


def save_gif(P, name='animation',
             T=None,
             fps=8,
             vmin=None, vmax=None,
             vwin=None,
             cmap='rainbow', dpi=150):
    '''
    Saves a [z,x,t] array as a gif.
    '''
    filename = name + '.gif'

    frames_div = 40

    nz, nx, nt = P.shape
    fig = plt.figure(dpi=dpi)
    fig.set_size_inches(2*nx/np.linalg.norm([nx, nz]),
                        2*nz/np.linalg.norm([nx, nz]))

    ax = fig.add_axes([0, 0, 1, 1], frameon=False, aspect=1)
    ax.set_xticks([])
    ax.set_yticks([])

    if vwin:
        for coord in ['xi', 'zi', 'xf', 'zf']:
            if not coord in vwin.keys():
                vwin[coord] = None
        vview = P[vwin['zi']:vwin['zf'], vwin['xi']:vwin['xf'], :]
        vmin = vview.min()
        vmax = vview.max()
    else:
        if not vmin:
            vmin = P.min()
        if not vmax:
            vmax = P.max()

    frames = []
    for f in range(1, nt, frames_div):
        frame = P[:, :, f]
        im = ax.imshow(frame, cmap, animated=True, vmin=vmin, vmax=vmax,
                       interpolation='nearest')
        frames.append([im])

    start = time()
    gif = manimation.ArtistAnimation(fig, frames, blit=True)
    gif.save(filename, writer='imagemagick', fps=fps)
    end = time()
    print('Gif saved; time spent:', end-start, end='.\n\n')
    plt.close(fig)
    return()


def save_gif(in_fn, out_fn):


    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    
    print(f'Reading array in "{filename}"')

    # opening read-only file as IO
    io = open(filename, "r")
    # read number of dimensions
    n = np.fromfile(io, dtype=np.int64, count=1)[0]
    # read dimensions as shape
    shape = np.fromfile(io, dtype=np.int64, count=n)
    # print(shape)
    # reading rest of file as matrix
    P = np.fromfile(io, dtype=np.float64,)
    # reshaping matrix as file
    P = P.reshape(shape, order='F')

    print('Reading array done!')

    print('saving snaps to file:')



if __name__ == '__main__':
    # defining parameters
    filename = "data/P.bin"

    print(f'Reading array in "{filename}"')

    # opening read-only file as IO
    io = open(filename, "r")
    # read number of dimensions
    n = np.fromfile(io, dtype=np.int64, count=1)[0]
    # read dimensions as shape
    shape = np.fromfile(io, dtype=np.int64, count=n)
    # print(shape)
    # reading rest of file as matrix
    P = np.fromfile(io, dtype=np.float64,)
    # reshaping matrix as file
    P = P.reshape(shape, order='F')

    print('Reading array done!')

    print('saving snaps to file:')

    plt.imshow(P[:, :, 201])

    save_gif(P, name='P_animation',
             fps=16,
             vmin=None, vmax=None,
             vwin=None,
             cmap='rainbow', dpi=150)
