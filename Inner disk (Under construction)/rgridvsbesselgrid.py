# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:51:30 2021

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


x = np.linspace(0,10,500)


plt.figure(figsize = (10,8))
for i in range(5):
    plt.plot(x, special.jv(i, x), label = '$J_%i$'%i)

plt.plot([0,10],[0,0], 'k-.')    
plt.legend(fontsize = 15)
plt.xlim(0,10)
plt.xlabel('$k$',fontsize = 15)
plt.ylabel('$J_q(k)$',fontsize = 15)