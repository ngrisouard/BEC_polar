# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:22:55 2020

@author: rundo
"""
'''
compute the first nzeros zeros of the bessel functions of the
norders first orders of the nkind-th kind and saves the result in the
current directory in dht.csv. Needs the routine bessel_zeros.m
'''
import csv
import numpy as np
from bessel_zeros import *
from scipy import special
import time

start = time.time()
norders = 1024
nzeros = 1024
nkind = 1

c = np.zeros((int(norders/2) + 1, nzeros+1), float)

c[0][:] = bessel_zeros(nkind, 0, nzeros+1)

for n in range(1, 100):
    print(bessel_zeros(nkind, n, nzeros)[0])
    print(n)
    
for n in range(100, 200):
    print(bessel_zeros(nkind, n, nzeros)[0])
    print(n)
    
for n in range(200, 300):
    print(bessel_zeros(nkind, n, nzeros)[0])
    print(n)

print(time.time() - start)
with open('dht.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    [writer.writerow(r) for r in c]