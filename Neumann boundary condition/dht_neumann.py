# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 19:35:01 2021

@author: rundo
"""

from numpy import genfromtxt
from scipy import special
import numpy as np

def dht_neumann(R, bessel_zeros, jmodes, ii, bessel_zeros_d, gain):
    
    C = bessel_zeros_d[abs(ii), jmodes-1]+gain[abs(ii)]
    #C = bessel_zeros[abs(ii), jmodes]
    c = bessel_zeros[abs(ii), :jmodes]
    
    Forward = np.zeros((jmodes, jmodes))
    Backward = np.zeros((jmodes, jmodes))
    
    if ii == 0:
        #c = np.insert(c,0,0)
        #c = c[:-1]
        for i in range(0, jmodes):
            for j in range(0, jmodes):
            
                Forward[j, i] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii, c[i])**2 
        for i in range(0, jmodes):
            for j in range(0, jmodes):
                Backward[i, j] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii, c[j])**2
    
    # elif ii == 1000:
    #     c = bessel_zeros[abs(ii), 1:jmodes+1]
    #     for i in range(0, jmodes):
    #         for j in range(0, jmodes):
            
    #             Forward[j, i] =  c[i] **2 *\
    #                 special.jv(ii, c[i]*c[j]/C) \
    #                     / (special.jv(ii, c[i])**2 * \
    #                        (c[i] **2- ii**2))
    
    
    #     for i in range(0, jmodes):
    #         for j in range(0, jmodes):
            
    #             Backward[i, j] =  c[j] **2 *\
    #                 special.jv(ii, c[i]*c[j]/C) \
    #                     / (special.jv(ii, c[j])**2 * \
    #                        (c[j] **2- ii**2))
    
    else:
        for i in range(1, jmodes):
            for j in range(1, jmodes):
            
                Forward[j, i] =  c[i] **2 *\
                    special.jv(ii, c[i]*c[j]/C) \
                        / (special.jv(ii, c[i])**2 * \
                           (c[i] **2- ii**2))


    
        for i in range(1, jmodes):
            for j in range(1, jmodes):
            
                Backward[i, j] =  c[j] **2 *\
                    special.jv(ii, c[i]*c[j]/C) \
                        / (special.jv(ii, c[j])**2 * \
                           (c[j] **2- ii**2))
    #Backward[0, 0] = 1
    #Backward[0,0] = 1/(4*np.pi*R**2/C**2)           
    r = R/C * c[:]
    k = c[:] / R
    Fi = 1/(np.pi * R**2)
    Bi = 4*np.pi*R**2/C**2
    '''   
    C = bessel_zeros_d[abs(ii), jmodes-1]
    c = bessel_zeros[abs(ii), :jmodes]
    
    r = R/C * c[:]
    k = c[:] / R
    
    I = abs(special.jv(ii, c))
    I[0] = 1/jmodes
    
    
    I = np.sqrt(2/C)/I
    I = np.outer(I, I) * special.jv(ii, np.outer(c, c/C))
    

    c[0] = 1/jmodes
    K = 2 * np.pi * R* abs(special.jv(ii, c)) * (1- ii**2/c[:]**2) /C   
    R = abs(special.jv(ii, c))* (1- ii**2/c[:]**2)/R
    
    K[0] = 1
    R[0] = 1
    '''
    #return k, r, I, K, R
    return r, k, Forward/(np.pi * R**2), Backward*4*np.pi*R**2/C**2#, Fi, Bi
    #return C, Forward, Backward
    
    
def dht_derichlet(R, bessel_zeros, jmodes, ii):
    
    C = bessel_zeros[abs(ii), jmodes]
    c = bessel_zeros[abs(ii), :jmodes]
    
    Forward = np.zeros((jmodes, jmodes))
    Backward = np.zeros((jmodes, jmodes))
    
    if ii == 0:
        for i in range(0, jmodes):
            for j in range(0, jmodes):
            
                Forward[j, i] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii+1, c[i])**2 

    
        for i in range(0, jmodes):
            for j in range(0, jmodes):
                Backward[i, j] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii+1, c[j])**2 
                    
    else:
        for i in range(1, jmodes):
            for j in range(1, jmodes):
            
                Forward[j, i] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii+1, c[i])**2 
        for i in range(1, jmodes):
            for j in range(1, jmodes):
                Backward[i, j] =  special.jv(ii, c[i]*c[j]/C) \
                    / special.jv(ii+1, c[j])**2
        Forward[0,0] = (np.pi * R**2)
        Backward[0,0] = 1/(4*np.pi*R**2/C**2)
    r = R/C * c[:]
    k = c[:] / R            
    #Backward[:, 0] = 1
    #Backward[:, 0] = 1
    return r, k, Forward/(np.pi * R**2), Backward*4*np.pi*R**2/C**2