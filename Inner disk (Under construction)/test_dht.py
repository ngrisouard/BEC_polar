# -*- coding: utf-8 -*-
"""
Created on Sat May 22 16:55:52 2021

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
from dht_annular import *
from dht import *
from ngobfft import *
from idht import *
import os

ii=100

jmodes = 128
bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')

R = 10

H,kk,rr,I,KK,RR, h = dht([],R,bessel_zeros,jmodes,ii)

I[0,0] = 1

hi = np.dot(I,I)


test = np.ones(jmodes)
#test = test*1j

test_k = np.dot(test / RR, I)*KK

test_r = np.dot(test_k / KK, I)*RR
