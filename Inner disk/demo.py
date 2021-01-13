# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 12:30:42 2021

@author: rundo
"""

import numpy as np
from scipy import special
from matplotlib import pyplot as plt
import csv
import time

def bessel(order, a, b, x):
    
    temp1 = special.jv(order, b * x) * special.yv(order, a * x)
    temp2 = special.jv(order, a * x) * special.yv(order, b * x)
    
    return temp1 - temp2


def bessel_zeros_ring(order, n_roots, a, b, accuracy):
    step = 0.1

    count = 0
    start = 0.5

    intervals = []

    while count <= n_roots:
        temp1 = bessel(order, a, b, start)
        temp2 = bessel(order, a, b, start+step)
    
        if np.sign(temp1) != np.sign(temp2):
            intervals.append([start, start+step])
            count += 1
        
        start += step

    roots = []  
    for r_range in intervals:
    
        left = np.float64(r_range[0])
        right = np.float64(r_range[1])
        mid = (left + right) / 2
    
        iteration = 0
        while abs(bessel(order, a, b, mid))/np.amax(abs(bessel(order, a, b, mid))) > accuracy: #and iteration <= 500:
        
            if np.sign(bessel(order, a, b, mid)) == np.sign(bessel(order, a, b, left)):
                left = np.copy(mid)
            else:
                right = np.copy(mid)
        
            mid = (left + right)/2
            iteration += 1
            if iteration%100 == 0:
                print(mid)
        #print(iteration)
        
        print(bessel(order, a, b, mid))
        roots.append(mid)
    
    return roots


order = 100
n_roots = 5
a = 1
b = 2
accuracy = 1e-14



roots = bessel_zeros_ring(order, n_roots, a, b, accuracy)
print('-------------------')
print(roots)

t = np.arange(roots[0]-1,roots[-1]+1,0.01)
f = bessel(order, a, b, t)
plt.plot(t, f)
