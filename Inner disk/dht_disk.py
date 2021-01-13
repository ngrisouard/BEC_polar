# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 19:51:44 2020

@author: rundo
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:17:12 2020

@author: rundo
%
%Quasi-discrete Hankel transform of integer order n.
%
%This implementation uses an array of Bessel roots, which
%is stored in 'dht.mat'. By default, the first 4097 roots
%of the Bessel functions of the first kind with order 0 to
%4 were precomputed. This allows DHT up to an order of 4
%with up to 4096 sampling points by default. Use JNROOTS
%for calculating more roots.
%
%Input:
% h      Function h(r)
%   or
% h      Signal h(r) *

% R      Maximum radius [m]
%   or
% R      Signal factors °                 {default}

%%%%%%%%%%%%%%%
c        Bessel zeros
%%%%%%%%%%%%%%%

% N      Number of sampling points        {256}
%   or
% K      Spectrum factors °               {default}


% n      Transform order                  {0}
%   or
% I      Integration kernel °             {default}
%
%Output:
% H      Spectrum H(k)
% k      Spatial frequencies [rad/m]
% r      Radial positions [m]
% I      Integration kernel
% K      Spectrum factors
% R      Signal factors
% h      Signal h(r)

% *)  Values request the presence of the kernel, samplings
%     and factors, which can be computed with empty input.
%
% °)  The computation of the integration kernel is slow com-
%     pared to the full transform. But if the factors and the
%     kernel are all present, they are reused as is.
"""
from numpy import genfromtxt
from scipy import special
import numpy as np
def bessel(order, a, b, x):
    
    temp1 = special.jv(order, b * x) * special.yv(order, a * x)
    temp2 = special.jv(order, a * x) * special.yv(order, b * x)
    
    return temp1 - temp2

def dht_disk(h, R0, R, c, N = 256, n = 0):
    k = None
    r = None
    if type(n) == list or type(n) == np.ndarray:
        if len(n) > 1:
            K = N
            I = n
            
    else:
        if (type(h) == list or type(h) == np.ndarray) and len(h) != 0:
            raise ValueError("need kernal")
            
        
        #c = genfromtxt('dht.csv', delimiter=',')
        #radius = int(R)
        #inner = int(N*R0/R)
        C = c[abs(n), N]
        c = c[abs(n), :N]
        
        #r_origin = R/C * c[:]
        r = c 
        
        k = c 
        
        I = np.sqrt(special.jv(n, c * R0)**2 - special.jv(n, c * R)**2) \
            / np.abs(c * special.jv(n, c*R))
        #if abs(n) > 0:
        #    I[0] = 1/N
        '''    
        
        '''
        K = I[:]
        R = I[:] 
        #I = np.sqrt(2) / I
        I = np.sqrt(2) / I
        #c_temp = np.copy(c)
        #cut_off = int(N*R0/radius)
        #c_temp[:cut_off] = 0
                #k[i] = 0
                
        apnapm = np.outer(c,c)
        apnb = c * R 
        
        jpyp1 = special.jv(n, apnapm)
        yp1 = special.yv(n, apnb)
        for i in range(N):
            jpyp1[i, :] = jpyp1[i, :] * yp1[i]
            
        jpyp2 = special.yv(n, apnapm)
        jp2 = special.jv(n, apnb)
        for i in range(N):
            jpyp2[i, :] = jpyp2[i, :] * jp2[i]
            
        
        
        I = np.outer(I, I) * (jpyp1 - jpyp2)
        #I[:cut_off, :] = 0
        #I[:, :cut_off] = 0
        '''
        if abs(n)>0:
            for i in range(N):
                if i < cut_off:
                    I[i][i] = 1
        '''
        
        
    if not any(h):
        H = h
    else:
        if type(h) != list and type(h) != float and type(h) != np.ndarray:
            h_result = h(r)
        else:
            h_result = h
        
        
        H = np.dot((h_result/ R), I) * K
        #H = np.dot(I, h_result/ R* K)
        
    return H, k, r, I, K, R, h
            