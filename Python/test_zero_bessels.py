# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 22:49:55 2020

@author: rundo
"""
import numpy as np
from bessel_zeros import *
from scipy import special
import mpmath as mp

jmodes = 0
nth = 5
R = 5

zero_bess_v = np.zeros((jmodes+1, nth))

#for i in range(0, jmodes+1):
for n in range(0, int(nth/2) + 1):
    zero_bess_v[:][n + int(nth/2)] = bessel_zeros(1, n, jmodes+1)

#for i in range(0, jmodes+1):
for n in range(1, int(nth/2)):
    zero_bess_v[:][n] = zero_bess_v[:][nth - n]
    
K_v = zero_bess_v/R

plottestv = np.zeros((jmodes+1,nth))

