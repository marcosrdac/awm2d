import numpy as np


def seisgain(d, dt=0.015, a=2., b=0.):
    nt, nx = d.shape
    t = np.arange(nt)*dt
    tgain = np.power(t, a) * np.exp(b*t)
    return tgain[:, None] * d.copy()
