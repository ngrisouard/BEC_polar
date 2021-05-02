# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 16:40:16 2021

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


#t = np.arange(0.5,10,0.01)
#f = bessel(128, 2, 10, t)
#plt.plot(t, f)

#freq = np.fft.fft(f)




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
        
        compare = r_range[0]-0.5
        left = np.float64(r_range[0])
        right = np.float64(r_range[1])
        mid = (left + right) / 2
    
        iteration = 0
        while abs(bessel(order, a, b, mid)/bessel(order, a, b, compare)) > accuracy and iteration <= 56:
        
            if np.sign(bessel(order, a, b, mid)) == np.sign(bessel(order, a, b, left)):
                left = np.copy(mid)
            else:
                right = np.copy(mid)
        
            mid = (left + right)/2
            iteration += 1
            #print(compare)
        
        print(bessel(order, a, b, mid))
        roots.append(mid)
    
    return roots


order = 128
n_roots = 256
accuracy = 1e-14

a = 8
b = 10
start = time.time()
c = []
for i in range(order+1):
    c.append(bessel_zeros_ring(i, n_roots, a, b, accuracy))
print(time.time() - start)
with open('dht_ring.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in c]