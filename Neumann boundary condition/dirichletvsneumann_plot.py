# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 21:52:31 2021

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
from dht_neumann import *
import os
import colorcet as cc



order = 2
r = np.linspace(0,20,1000)

f = special.jv(order, r)


plt.figure(figsize = (10,7))
plt.plot(r,f,'k', linewidth = 1)
plt.plot(r,np.zeros(1000),'-.k', linewidth = 1)


n_root = np.genfromtxt('dht_neumann.csv', delimiter=',')
d_root = np.genfromtxt('dht.csv', delimiter=',')


#plt.scatter(n_root[order, 4], special.jv(order, n_root[order, 4]), c='r', marker='x')

plt.scatter(d_root[order, :6], special.jv(order, d_root[order, :6]),s=100, c='b', marker='s', label='$k^{Dir}$, Dirichlet roots')
plt.scatter(n_root[order, :6], special.jv(order, n_root[order, :6]),s=100, c='r', marker='^', label='$k^{Neu}$, Neumann roots')

plt.plot([n_root[order, 5], n_root[order, 5]], [special.jv(order, n_root[order, 5]), 0.3], 'b', linewidth = 0.5)
plt.plot([n_root[order, 4], n_root[order, 4]], [special.jv(order, n_root[order, 4]), 0.3], 'b', linewidth = 0.5)
plt.plot([d_root[order, 4], d_root[order, 4]], [special.jv(order, d_root[order, 4]), 0.4], 'b', linewidth = 0.5)
plt.plot([d_root[order, 5], d_root[order, 5]], [special.jv(order, d_root[order, 5]), 0.4], 'b', linewidth = 0.5)
plt.plot([d_root[order, 0], d_root[order, 0]], [special.jv(order, d_root[order, 0]), -0.3], 'b', linewidth = 0.5)

plt.annotate('$k_{q,N+1}^{Neu}$', [n_root[order, 5], 0.32], fontsize = 18)
plt.annotate('$k_{q,N}^{Neu}$', [n_root[order, 4], 0.32], fontsize = 18)
plt.annotate('$k_{q,N}^{Dir}$', [d_root[order, 4], 0.42], fontsize = 18)
plt.annotate('$k_{q,N+1}^{Dir}$', [d_root[order, 5], 0.42], fontsize = 18)
plt.annotate('$k_{q,1}^{Dir}$ & $k_{q,1}^{Neu}$', [d_root[order, 0], -0.32], fontsize = 18)
plt.xlabel('$k$')
plt.ylabel('$J_q(k)$')
plt.legend(fontsize = 13, loc = 4)