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


def dht(h, R, c, N = 256, n = 0):
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
            