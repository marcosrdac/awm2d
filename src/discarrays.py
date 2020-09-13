import numpy as np

def discarray(filename, mode='r', dtype=float, shape=None, order='C'):
    file_mode = f'{mode[0]}b{mode[1:]}'
    if not isinstance(shape, tuple):
        shape = (shape,)
    with open(filename, file_mode) as io:
        if 'w' in mode:
            ndims_shape = np.array((len(shape), *shape), dtype=np.int64)
            ndims_shape.tofile(io)
        if 'r' in mode:
            ndims = np.fromfile(io, dtype=np.int64, count=1)[0]
            shape = tuple(np.fromfile(io, dtype=np.int64, count=ndims))
        offset = io.tell()
        arr = np.memmap(io, dtype=dtype, mode=mode,
                        offset=offset, shape=shape, order=order)
    return(arr)
