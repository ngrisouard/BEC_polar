# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 21:52:44 2020

@author: rundo
"""
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from scipy import interpolate
from dht_disk import *
from dht_annular import *
from dht import *
from ngobfft import *
from idht import *
import os

def bessel(order, R0, R, x, root):
    
    temp1 = special.jv(order, root * x /R) * special.yv(order, root*R0/R)
    temp2 = special.jv(order, root*R0/R ) * special.yv(order, root * x /R)
    
    return temp1 - temp2

#%%
ii=12

jmodes = 128
R0 = 1
R = 10
nr = 256
r = np.linspace(R0, R, nr)

bessel_zeros = np.genfromtxt('dht_ring_10_1.csv', delimiter=',')

Forward, Backward, KK, RR = dht_annular(r, R0, R, bessel_zeros, jmodes, ii)
identity = np.dot(Forward, Backward)* np.pi**2 / 2
#%%
#order = 40
test = 1*(bessel(ii, R0, R, RR, bessel_zeros[ii,1]) + 1*bessel(ii, R0, R, RR, bessel_zeros[ii,3]) + 1.5*bessel(ii, R0, R, RR, bessel_zeros[ii,10])) #* np.exp(0j*np.pi/2)
#test = np.ones(nr)
test[0] = 0
test[-1] = 0
k = np.dot(test, Forward)
#%%
dt=0.005
#k = np.dot(f, Forward)

#f = np.dot(k*np.exp(-1j*0.5*(KK**2)*dt)* np.pi**2 / 2, Backward)
f = np.dot(k*np.pi**2 / 2, Backward)
k = np.dot(f, Forward)


plt.plot(RR, f, label='transformed')
plt.plot(RR, test, label='input')
plt.title('Order = %i'%ii)
plt.legend()
print('error:',np.sum(np.abs(f - test)))

#%%

ii=0

jmodes = 256
bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')

R = 10

jmodes = 256


r = np.linspace(0, R, jmodes)


dx = R/jmodes


Forward_disk = np.zeros((jmodes, jmodes))

for i in range(jmodes):
    for j in range(jmodes):
        #k[i] += r[j]*dx*bessel(ii, R, r[j], bessel_zeros[ii, i])*test[j]
        
        Forward_disk[j,i] = r[j]*dx*special.jv(ii, r[j]*bessel_zeros[ii, i]/R)
        


Backward_disk = np.zeros((jmodes, jmodes))


for i in range(jmodes):
    for j in range(jmodes):
        #f[i] += (bessel_zeros[ii, j]/R)**2 * special.jv(ii, bessel_zeros[ii, j])**2 \
        #    * k[j] * bessel(ii, R, r[i], bessel_zeros[ii, j]) \
        #        / (special.jv(ii, bessel_zeros[ii, j]*R0/R)**2 - special.jv(ii, bessel_zeros[ii, j])**2)
        if j == 0:
            Backward_disk[j,i] = 0
        else:
            Backward_disk[j,i] = special.jv(ii, r[i]*bessel_zeros[ii, j]/R) / special.jv(ii+1, bessel_zeros[ii, j])**2
        
identity_disk = np.dot(Forward_disk, Backward_disk)




test = special.jv(ii, r*bessel_zeros[ii, 2]/R) #bessel(ii, R, r, bessel_zeros[ii,1]) + 10*bessel(ii, R, r, bessel_zeros[ii,3]) + bessel(ii, R, r, bessel_zeros[ii,6])
#test = np.ones(jmodes)
#test[0] = 0
#test[-1] = 0


k = np.dot(test, Forward_disk)

f = np.dot(k*(2/R**2), Backward_disk)
                

plt.plot(r, f)
plt.plot(r, test)


#%%

ii=0

jmodes = 256
bessel_zeros = np.genfromtxt('dht_ring.csv', delimiter=',')


#k, r, iI_f, iI_b, kK, rR, J_q, Y_q, scale = dht_annular([], 1, 1.5, bessel_zeros, jmodes, ii)

#hi = np.dot(iI_f, iI_b)


R0 = 1
R = 1.25

r = np.linspace(R0, R, jmodes)


dx = (R-R0)/jmodes



Forward = np.zeros((jmodes, jmodes))

for i in range(jmodes):
    for j in range(jmodes):
        #k[i] += r[j]*dx*bessel(ii, R, r[j], bessel_zeros[ii, i])*test[j]
        
        Forward[j,i] = r[j]*dx*bessel(ii, R, r[j], bessel_zeros[ii, i])
        


Backward = np.zeros((jmodes, jmodes))


for i in range(jmodes):
    for j in range(jmodes):
        #f[i] += (bessel_zeros[ii, j]/R)**2 * special.jv(ii, bessel_zeros[ii, j])**2 \
        #    * k[j] * bessel(ii, R, r[i], bessel_zeros[ii, j]) \
        #        / (special.jv(ii, bessel_zeros[ii, j]*R0/R)**2 - special.jv(ii, bessel_zeros[ii, j])**2)
                
        Backward[j,i] = (bessel_zeros[ii, j]/np.sqrt(R))**2 * special.jv(ii, bessel_zeros[ii, j])**2 \
             * bessel(ii, R, r[i], bessel_zeros[ii, j]) \
                / (special.jv(ii, bessel_zeros[ii, j]*R0/R)**2 - special.jv(ii, bessel_zeros[ii, j])**2)

        #Backward[j,i] = 2*bessel(ii, R, r[i], bessel_zeros[ii, j])/(R**2 *(bessel(ii, R, R, bessel_zeros[ii, j])**2 -\
        #                                                                   bessel(ii+1, R, R, bessel_zeros[ii, j])*bessel(ii-1, R, R, bessel_zeros[ii, j]))\
        #                 - R0**2 * (bessel(ii, R, R0, bessel_zeros[ii, j])**2 - bessel(ii+1, R, R0, bessel_zeros[ii, j])*bessel(ii-1, R, R0, bessel_zeros[ii, j])))

identity = np.dot(Forward, Backward)* np.pi**2 / 2








