# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:51:54 2021

@author: rundo
"""

from numpy import genfromtxt
from scipy import special
import numpy as np

def bessel(order, R0, R, x, root):
    
    temp1 = special.jv(order, root * x /R) * special.yv(order, root*R0/R)
    temp2 = special.jv(order, root*R0/R ) * special.yv(order, root * x /R)
    
    return temp1 - temp2

def dht_annular(r, R0, R, bessel_zeros, jmodes, ii):
    
    #r = np.linspace(R0, R, jmodes)


    dx= r[1]-r[0]



    Forward = np.zeros((len(r), jmodes))

    for i in range(jmodes):
        for j in range(len(r)):
        #k[i] += r[j]*dx*bessel(ii, R, r[j], bessel_zeros[ii, i])*test[j]
        
            Forward[j,i] = r[j]*dx*bessel(ii,R0, R, r[j], bessel_zeros[abs(ii), i])


    Backward = np.zeros((jmodes, len(r)))


    for i in range(len(r)):
        for j in range(jmodes):
        #f[i] += (bessel_zeros[ii, j]/R)**2 * special.jv(ii, bessel_zeros[ii, j])**2 \
        #    * k[j] * bessel(ii, R, r[i], bessel_zeros[ii, j]) \
        #        / (special.jv(ii, bessel_zeros[ii, j]*R0/R)**2 - special.jv(ii, bessel_zeros[ii, j])**2)
                
            Backward[j,i] = (bessel_zeros[abs(ii), j]/R)**2 * special.jv(ii, bessel_zeros[abs(ii), j])**2 \
             * bessel(ii, R0, R, r[i], bessel_zeros[abs(ii), j]) \
                / (special.jv(ii, bessel_zeros[abs(ii), j]*R0/R)**2 - special.jv(ii, bessel_zeros[abs(ii), j])**2)
            
            #Backward[j,i] = -2*bessel(ii, R, r[i], bessel_zeros[abs(ii), j]) \
            #                / ((bessel_zeros[abs(ii), j])**2 * (R**2 * (bessel(ii, R, R, bessel_zeros[abs(ii), j])**2 - bessel(ii-1, R, R, bessel_zeros[abs(ii), j])*bessel(ii+1, R, R, bessel_zeros[abs(ii), j]))\
            #                   - R0**2 * (bessel(ii, R, R0, bessel_zeros[abs(ii), j])**2 - bessel(ii-1, R, R0, bessel_zeros[abs(ii), j])*bessel(ii+1, R, R0, bessel_zeros[abs(ii), j]))))
            
            
            
    return Forward, Backward, bessel_zeros[ii, :jmodes]/R, r
    