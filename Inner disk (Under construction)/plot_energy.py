# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:19:27 2021

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

new = io.loadmat('./output/kernal/energy_new')
old = io.loadmat('./output/kernal/energy_old')

newL = io.loadmat('./output/kernal/L_new')
oldL = io.loadmat('./output/kernal/L_old')

enew = new['E']
eold = old['E']

lnew = newL['L']
lold = oldL['L']

start = 0
end = 501

it = range(len(enew[0]))

plt.figure(figsize = (8,8))
plt.plot(it, enew[0], label = 'New method')
plt.plot(it, eold[0], label = 'Old method')
plt.xlabel('Iterations', fontsize = 20)
plt.ylabel('E', fontsize = 20)
plt.grid('on')
plt.xlim(0,500)
plt.ylim(-45.2,-44.8)
plt.legend()



plt.figure(figsize = (8,8))
plt.plot(it, np.imag(lnew[0]), label = 'New method')
plt.plot(it, np.imag(lold[0]), label = 'Old method')
plt.xlabel('Iterations', fontsize = 20)
plt.ylabel('L', fontsize = 20)
plt.grid('on')
plt.xlim(0,300)
plt.legend()