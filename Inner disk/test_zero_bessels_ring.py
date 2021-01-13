# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 15:04:17 2021

@author: rundo
"""

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from scipy import interpolate
from dht_disk import *
from dht import *
from ngobfft import *
from idht import *
import os
def bessel(order, a, b, x):
    
    temp1 = special.jv(order, b * x) * special.yv(order, a * x)
    temp2 = special.jv(order, a * x) * special.yv(order, b * x)
    
    return temp1 - temp2

bessel_zeros = np.genfromtxt('dht_ring.csv', delimiter=',')
alala = np.genfromtxt('dht.csv', delimiter=',')

a = 8
b = 10
zeros = np.zeros(np.shape(bessel_zeros))
for i in range(np.shape(bessel_zeros)[0]):
    zeros[i, :] = bessel(i, a, b, bessel_zeros[i, :])    
