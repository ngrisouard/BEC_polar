# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 22:49:55 2020

@author: rundo
"""
import numpy as np
from bessel_zeros import *
from scipy import special
import mpmath as mp
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

jmodes = 30
nth = 30
R = 5

zero_bess_v = np.zeros((jmodes+1, nth))

for i in range(0, jmodes+1):
    for n in range(0, int(nth/2)+1):
        zero_bess_v[i][n + int(nth/2) - 1] = bessel_zeros(1, n, jmodes+1)[i]

for i in range(0, jmodes+1):
    for n in range(0, int(nth/2)-1):
        zero_bess_v[i][n] = zero_bess_v[i][nth - n - 2]

K_v = zero_bess_v/R

plottestv = np.zeros((jmodes+1,nth))

for i in range(0, jmodes+1):
    for n in range(-int(nth/2) + 1, int(nth/2)+1):
        plottestv[i][n-1 + int(nth/2)] = abs(special.jv(n, zero_bess_v[i][n-1 + int(nth/2)]))
        
        

plt.pcolormesh(plottestv, norm=LogNorm(vmin=1e-18, vmax=plottestv.max()))
plt.colorbar()
