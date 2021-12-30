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
    
    temp1 = special.jv(order, x * a/b) * special.yv(order,  x)
    temp2 = special.jv(order, x) * special.yv(order, x * a/b)
    
    return temp1 - temp2


def swap(R0, R, k, zeros, order):
    
    return (R*(k-zeros[order,0]) - R0 * (k-zeros[order,-1])) / (zeros[order,-1]-zeros[order,0])

def swap_2(R0, R, k, zeros, order):
    
    return k/(zeros[order,-1]-zeros[order,0])*(R-R0)

bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')
ring_zeros = np.genfromtxt('dht_ring.csv', delimiter=',')

R0 = 8
R = 10
order = 1
zeros = np.zeros(np.shape(ring_zeros))


root_grid = swap_2(R0, R, ring_zeros[order,:], ring_zeros, order)


