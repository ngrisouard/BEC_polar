# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 18:14:29 2020

@author: rundo
"""

'''
compare matlab zeros with python zeros
'''
import scipy.io
from numpy import genfromtxt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

python = genfromtxt('dht.csv', delimiter=',')
matlab = scipy.io.loadmat('dht.mat')
mat = matlab['c']

diff = np.zeros((513, 1025))

for i in range(513):
    for j in range(1025):
        diff[i][j] = abs(python[i][j] - mat[i][j])

plt.pcolormesh(diff, norm=LogNorm(vmin=1e-18, vmax=diff.max()))
plt.colorbar()