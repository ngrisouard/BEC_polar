# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 01:22:25 2021

@author: rundo
"""

from matplotlib import pyplot as plt
from scipy import special
import numpy as np

def bessel(order, R0, R, x, root):
    
    temp1 = special.jv(order, root * x /R0) * special.yv(order, root * R/R0)
    temp2 = special.jv(order, root * R/R0) * special.yv(order, root * x /R0)
    
    return temp1 - temp2


def bessel2(order, R0, R, x, root):
    
    temp1 = special.jv(order, root * x) * special.yv(order, root * R)
    temp2 = special.jv(order, root * R) * special.yv(order, root * x)
    
    return temp1 - temp2

#bessel_zeros = np.genfromtxt('dht.csv', delimiter=',')
ring_zeros = np.genfromtxt('dht_ring_10_1.csv', delimiter=',')


#R = 10
#r = np.linspace(0,R, 100)
typelist = ['k-','k--','k:','k-.']
#order = [0,1,4,6]
'''
plt.figure(figsize = (8,8))



for i in range(1,2):
    for j in order:
        lab = 'J$_{%i}$'%j + '(k$_{%i,}$'%j+'$_{%i}r/R$)'%i
        jj = special.jv(j,r*bessel_zeros[j,i]/R)

        plt.plot(r/R,jj, typelist[order.index(j)], label=lab)
        #plt.plot(r/R,jj, typelist[i-1], label=lab)

plt.plot(r/R, np.zeros(len(r)), 'k-.', linewidth=1)
plt.xlim(0,1)
plt.xlabel('r/R', fontsize=15)
plt.legend(fontsize=15)
plt.show()
'''
R0 = 1
R = 10

r = np.linspace(R0,R, 500)
order = 1


plt.figure(figsize = (8,8))
rootlist = [0,2,3,4]

for i in rootlist:
    lab = 'k$_{%i,}$'%order+'$_{%i}$'%(i+1)
    plt.plot(r, bessel2(order, R0,R, r, ring_zeros[abs(order),i]), typelist[rootlist.index(i)], label=lab)
#plt.plot(r, bessel(order, R, r, ring_zeros[order,2]), 'k--')
#plt.plot(r, bessel(order, R, r, ring_zeros[order,4]), 'k:')
#plt.plot(r, bessel(order, R, r, ring_zeros[order,5]), 'k-.')
plt.grid('on')
plt.legend(fontsize=15)
plt.title('Order: %i'%order, fontsize=15)
#plt.ylim(-0.1,0.1)
#plt.xlim(R0,R)
plt.xlabel('r', fontsize=15)
plt.show()


'''
order = [0,1,4,6]

plt.figure(figsize = (8,8))
r = np.linspace(R0,R, 100)
typelist = ['k-','k--','k:','k-.']

for i in range(1,2):
    for j in order:
        lab = 'J$_{%i}$'%j + '(k$_{%i,}$'%j+'$_{%i}r/R$)'%i
        jj = bessel(j,R , r, ring_zeros[j,i])
        plt.plot(r,jj, typelist[order.index(j)], label=lab)
        #plt.plot(r/R,jj, typelist[i-1], label=lab)
        
#plt.plot(r, np.zeros(len(r)), 'k-.', linewidth=1)
plt.xlim(8,10)
plt.xlabel('r', fontsize=15)
plt.legend(fontsize=15)
plt.show()

x = np.linspace(1, 15, 1000)
plt.figure(figsize = (12,8))
plt.plot(x, rootfunc(1, R0, R, x), 'k-')
plt.grid('on')
plt.xlabel('k', fontsize = 15)
plt.show()
'''