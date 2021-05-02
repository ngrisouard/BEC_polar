# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 21:52:44 2020

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

ii=100
ker = io.loadmat('output/kernal/ker_%i.mat'%(ii))
I = ker['I']
kk = ker['kk']
KKa = ker['KK']
rr = ker['rr']
RRa = ker['RR']
bessel_zeros = np.genfromtxt('dht_ring.csv', delimiter=',')