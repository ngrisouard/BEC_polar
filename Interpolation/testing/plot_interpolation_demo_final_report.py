# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 20:43:41 2021

@author: rundo
"""

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from scipy import interpolate
from dht import *
from ngobfft import *
from idht import *
from interpolation import *
from conservative import *
import os


R = 10


nr = 11        #radial points
dr = R / (nr - 1)
r = np.arange(0, R+dr, dr)


bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')


q = 0
roots = bessel_zeros[1, :11]
maxroot = bessel_zeros[1, 11]

root_grid = roots*R/maxroot


y = np.sin(root_grid)*np.exp(0.2*root_grid)
y2 = np.sin(r)*np.exp(0.2*r)
plt.figure(figsize = (16,8))
plt.scatter(root_grid, y, color='r', label = 'Bessel roots grid, $k_{q,j}R/k_{q,N+1}$', marker='x', s=100)
plt.scatter(r,y2, label = 'Evenly sampled $r$ grid',s=60)

x = np.linspace(0,10,200)
line = np.sin(x)*np.exp(0.2*x)

plt.plot(x, line, 'k', linewidth = 1, label='Signal $f(r)$')
plt.plot([0,10],[0,0],'k-.')

for i in range(len(root_grid)):
    plt.plot([root_grid[i],root_grid[i]],[0, y[i]],'r', linewidth=1)
    plt.plot([r[i],r[i]],[0, y2[i]],'b', linewidth=1)


plt.xlim(-0.1,10.1)
plt.legend(fontsize = 20)
plt.yticks(np.arange(1), " ")
