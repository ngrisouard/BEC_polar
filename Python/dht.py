# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:17:12 2020

@author: rundo
"""
from numpy import genfromtxt
from scipy import special
import numpy as np


def dht(h, R, N = 256, n = 0):
    if type(n) == list or type(n) == np.ndarray:
        if len(n) > 1:
            K = N
            I = n
            
    else:
        if type(h) == list and len(h) != 0:
            raise ValueError("need kernal")
            
        
        c = genfromtxt('dht.csv', delimiter=',')
        C = c[n, N]
        c = c[n, :N]
        r = R/C * c[:]
        k = c[:] / R
        I = abs(special.jv(1+n, c))
        if n > 0:
            I[0] = 1/N
            
        K = 2 * np.pi * R/C * I[:]
        R = I[:]/R
        I = np.sqrt(2/C) / I
        I = np.outer(I, I) * special.jv(n, np.outer(c/C, c))
        
    if not h:
        H = h
    else:
        if type(h) != list and type(h) != float:
            h_result = h(r)
        else:
            h_result = h
        
        H = I*(h_result/ R) * K
        
    return H, k, r, I, K, R, h
            