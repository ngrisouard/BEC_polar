# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 16:45:28 2020

@author: rundo
"""

from dht import dht
from scipy import special
from numpy import genfromtxt
from scipy import io
import numpy as np
import csv
from dht_neumann import *
from matplotlib import pyplot as plt


def bessr(d, a, x):
    if d == 1:
        return special.jv(a, x) / special.jv(a+1, x)
    elif d == 2:
        return special.yv(a, x) / special.yv(a+1, x)
    elif d == 3:
        return a/x - special.jv(a+1, x) / special.jv(a, x)
    else:
        return a/x - special.yv(a+1, x) / special.yv(a, x)
#c = genfromtxt('dht.csv', delimiter=',')

R = 10
jmodes = 128

ii = 10
N=256
h=[]
dt = 5e-2

bessel_zeros_neu = np.genfromtxt('dht_neumann.csv', delimiter=',')

bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')

H,kk,rr,I,KK,RR, h = dht([],R,bessel_zeros,jmodes,ii)
if ii != 0:
    I[0,0]=1

k_n, rr, Forward, Backward = dht_neumann(R, bessel_zeros_neu, jmodes, ii, bessel_zeros)

r_d,k_d, Forward_d, Backward_d = dht_derichlet(R, bessel_zeros, jmodes, ii)
#save_dict = {'kk':kk, 'rr':rr, 'KK':KK, 'RR':RR, 'I':I}

#io.savemat('./output/ker.mat', save_dict)

#read = io.loadmat('./output/ker_3.mat')

#I = read['I']
result = np.dot(Forward, Backward)


result2 = np.dot(Forward_d, Backward_d)

#result3 = np.dot(I,I)

root_grid = bessel_zeros[ii, :jmodes]*R\
                    /bessel_zeros[ii, jmodes]

#f = np.tanh((R-root_grid)/np.sqrt(2)) #* np.exp(0j*Thet)#special.jv(ii, bessel_zeros[ii, 3]*root_grid)#np.ones(jmodes)

f = np.ones(jmodes)

k = np.dot(Backward, f) #- f[0]/(np.pi * R**2)

f_res = np.dot(Forward, k*np.exp(-0j*0.5*(k_n**2)*dt))

k2 = np.dot(f/RR, I)*KK
f_res2 = np.dot(k2*np.exp(-1j*0.5*(k_d**2)*dt)/KK, I)*RR


plt.plot(f**2)
plt.plot(np.abs(f_res)**2)
#plt.plot(np.abs(f_res2)**2)

#sub = bessr(3, n, bessel_zeros_neu[n,:])
'''
#write
with open('./output/ker.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(save_dict.keys())
    writer.writerows(zip(*save_dict.values()))
        
#read

input_file = csv.DictReader(open("./output/ker.csv"))
for row in input_file:
    print(row)

'''

    
