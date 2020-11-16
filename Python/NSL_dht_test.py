# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:02:58 2020

@author: rundo
"""

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from dht import *
from ngobfft import *
import os
from scipy import interpolate

#ker = io.loadmat('output/kernal/ker_0.mat')
#misc = io.loadmat('output/misc.mat')
#wf = io.loadmat('output/psi_00000000.mat')

#kt = k_of_x(Thet[0]) 
#fr = obfft(Thet[0], wf, -1) #fft to isolate every angular mode


#r2 = k_of_x(kt)
'''
x = [1,2,3,4]
y = [3,4,6,1]
sb = interpolate.interp1d(x, y,kind = 'cubic')
print(sb(1.5))
'''

x = np.array([[1,2,3,4],
              [4,3,2,1],
              [5,6,7,8],
              [4,5,6,7]])

y = np.array([[2,2,2,2],
              [2,2,2,2],
              [2,2,2,2],
              [2,2,2,2]])
