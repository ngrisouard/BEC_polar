# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 21:04:51 2020

@author: rundo
"""
from scipy import special
import mpmath as mp


def bessel_zeros(d, a, n):
    mp.dps = 15
    if d == 1:      #Ja
        zeros = []
        for i in range(1, n+1):
            zeros.append(float(mp.besseljzero(a, i)))
        
        return zeros
        
    elif d == 2:    #Ya
        zeros = []
        for i in range(1, n+1):
            zeros.append(float(mp.besselyzero(a, i)))
        
        return zeros
    elif d == 3:    #Ja'
        zeros = []
        for i in range(1, n+1):
            zeros.append(float(mp.besseljzero(a, i, derivative = 1)))
        
        return zeros
    else:    #Ya'
        zeros = []
        for i in range(1, n+1):
            zeros.append(float(mp.besselyzero(a, i, derivative = 1)))
        
        return zeros