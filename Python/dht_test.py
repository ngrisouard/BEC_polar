# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 16:45:28 2020

@author: rundo
"""

from dht import dht
from scipy import special
from numpy import genfromtxt
import numpy as np

#c = genfromtxt('dht.csv', delimiter=',')

R = 10
jmodes = 128
n = 215
ii = n
N=256
h=[]

'''
if type(n) == list:
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
    if type(h) != list or type(h) != float:
        h_result = h(r)
    else:
        h_result = h
        
    H = I*(h_result/ R) * K
'''




H, k, r, I, K, R, h = dht([],R,jmodes,ii)