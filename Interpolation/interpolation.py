# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 17:56:31 2021

@author: rundo


r: radial grids
R: max radius
n: transform order
bessel_zeros: precalculated bessel roots


"""





from numpy import genfromtxt
from scipy import special
import numpy as np

def compute_interp(r, R, n, bessel_zeros, N):
    
    jnk = bessel_zeros[abs(n), :N]
    
    W = bessel_zeros[abs(n), N]
    
    
    I = np.ndarray((len(jnk),(len(r))))
    

    for j in range(len(r)):
            
        denominator = (special.jn(n+1, jnk)*(jnk**2-(r[j]*W/R)**2))
        if abs(n) >= 0:
            denominator[0] = 1/N
            
        I[:,j] = 2*jnk*special.jn(n, r[j]*W/R)/denominator
    
    
    
    return I




def compute_interp_alt(r, R, n, bessel_zeros, N):
    
    if n != 0:
        jnk = bessel_zeros[abs(n), 1:N+1]
    
        W = bessel_zeros[abs(n), N+1]
    else:
        jnk = bessel_zeros[abs(n), :N]
    
        W = bessel_zeros[abs(n), N]
    
    
    I = np.ndarray((len(jnk),(len(r))))
    

    for j in range(len(r)):
            
        denominator = (special.jn(n+1, jnk)*(jnk**2-(r[j]*W/R)**2))
        #if abs(n) > 0:
        #    denominator[0] = 1/N
            
        I[:,j] = 2*jnk*special.jn(n, r[j]*W/R)/denominator
    
    
    return I


