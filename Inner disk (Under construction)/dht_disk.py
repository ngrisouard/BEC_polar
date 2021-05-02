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
            
        
        #radius = int(R)
        #inner = int(N*R0/R)
        C = c[abs(n), N]
        c = np.array(c[abs(n), :N],complex)
        
        r = R/C * c[:]
        k = c[:] / R
        
        scale = np.abs(c[:]) * np.abs(special.jv(n, c)) \
            / (np.sqrt(special.jv(n, c*R0/R)**2 - special.jv(n, c)**2) *R)
        
        K = 2  * R  / (C * scale) * np.pi
        R = 1 / (R * scale)
        
        I = np.sqrt(2/C) * scale
        I = np.outer(I, I)
        
        generalJ = special.jv(n, np.outer(c, c/C))
        generalY = special.yv(n, np.outer(c, c/C))
        for i in range(N):
            generalJ[:, i] = generalJ[:, i] * special.yv(n, c)
            generalY[:, i] = generalY[:, i] * special.jv(n, c)
        
        I = I * (generalJ - generalY)
        
        
        
        
        
    if not any(h):
        H = h
    else:
        if type(h) != list and type(h) != float and type(h) != np.ndarray:
            h_result = h(r)
        else:
            h_result = h
        
        
        H = np.dot(h_result, I) * K
        #H = np.dot(I, h_result/ R* K)
        
    return H, k, np.abs(r), I, K, R, h
            