# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 09:49:27 2021

@author: rundo
"""

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from scipy import interpolate
from scipy import special
#from dht import *
#from ngobfft import *
#from idht import *
#from interpolation import *
#from conservative import *
#from dht_neumann import *
import os
import colorcet as cc




def Forward_D(q, R, k, N):
    S = k[abs(q), N]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = 1/special.jv(q+1, k[j])**2 \
                    * special.jv(q, k[j]*k[i]/S)
    
    else:
        for i in range(1,N):
            for j in range(1,N):
                matrix[i,j] = 1/special.jv(q+1, k[j])**2 \
                    * special.jv(q, k[j]*k[i]/S)
        #matrix[0,0] = 1
    return 4*np.pi*R**2 / S**2 * matrix


def Backward_D(q, R, k, N):
    S = k[abs(q), N]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = 1/special.jv(q+1, k[j])**2 \
                    * special.jv(q, k[j]*k[i]/S)
    
    else:
        for i in range(1,N):
            for j in range(1,N):
                matrix[i,j] = 1/special.jv(q+1, k[j])**2 \
                    * special.jv(q, k[j]*k[i]/S)
        #matrix[0,0] = 1
    return 1/(np.pi * R**2) * matrix

def Derivative_D(q, R, k, N):#, neu):
    #S = neu[abs(q), N]
    S = k[abs(q), N]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = k[j] / special.jv(q+1, k[j])**2 \
                    * (special.jv(q-1, k[j]*k[i]/S) - special.jv(q+1, k[j]*k[i]/S))
    else:
        for i in range(1,N):
            for j in range(1,N):
                matrix[i,j] = k[j] / special.jv(q+1, k[j])**2 \
                    * (special.jv(q-1, k[j]*k[i]/S) - special.jv(q+1, k[j]*k[i]/S))
                    
    return 1/(2*np.pi * R**3) * matrix

def multiR_D(q, R, k, N):
    S = k[abs(q), N]
    matrix = np.zeros((N,N))
    for i in range(N):
        matrix[i,i] = k[abs(q),i]*R/S
    
    return matrix


def divR_D(q, R, k, N):
    '''
    S = k[abs(q), N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = k[abs(q),j] / special.jv(q+1, k[abs(q),j])**2 \
                    * (special.jv(q-1, k[abs(q),j]*k[abs(q),i]/S) + special.jv(q+1, k[abs(q),j]*k[abs(q),i]/S)) /q
    else:
        for i in range(1,N):
            for j in range(1,N):
                matrix[i,j] = k[abs(q),j] / special.jv(q+1, k[abs(q),j])**2 \
                    * (special.jv(q-1, k[abs(q),j]*k[abs(q),i]/S) + special.jv(q+1, k[abs(q),j]*k[abs(q),i]/S)) /q
    
    
    return 1/(2*np.pi * R**3) * matrix
    '''
    S = k[abs(q), N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            matrix[i,i] = S/(k[abs(q),i]*R)
    
    else:
        for i in range(1,N):
            matrix[i,i] = S/(k[abs(q),i]*R)
    
        #matrix[0,0] = 1/N
    
    return matrix
    
def dTheta_D(q, R, k, N):
    matrix = np.zeros((N,N),complex)
    for i in range(N):
        matrix[i,i] = 1j*q
    
    return matrix


def Laplacian_D(q, R, k, N):
    matrix = np.zeros((N,N))
    for i in range(N):
        matrix[i,i] = -(k[abs(q),i]/R)**2
    return matrix

def Forward_N(q, R, k, N, kd):
    S = kd[abs(q), N-1]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = 1 / special.jv(q, k[j])**2  \
                    * special.jv(q, k[j]*k[i]/S)
                    
    else:
        for i in range(1, N):
            for j in range(1,N):
                matrix[i,j] = k[j]**2 / ((k[j]**2 - q**2)* special.jv(q, k[j])**2) \
                    * special.jv(q, k[j]*k[i]/S)
                    
    return 4*np.pi*R**2 / S **2 * matrix

def Backward_N(q, R, k, N, kd):
    S = kd[abs(q), N-1]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = 1 / special.jv(q, k[j])**2  \
                    * special.jv(q, k[j]*k[i]/S)
                    
    else:
        for i in range(1, N):
            for j in range(1,N):
                matrix[i,j] = k[j]**2 / ((k[j]**2 - q**2)* special.jv(q, k[j])**2) \
                    * special.jv(q, k[j]*k[i]/S)
                    
    return 1 / (np.pi * R**2) * matrix


def Derivative_N(q, R, k, N, kd):#, neu):
    #S = neu[abs(q), N]
    S = kd[abs(q), N-1]
    k = k[abs(q), :N]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            for j in range(N):
                matrix[i,j] = k[j] / special.jv(q, k[j])**2 \
                    * (special.jv(q-1, k[j]*k[i]/S) - special.jv(q+1, k[j]*k[i]/S))
    else:
        for i in range(1,N):
            for j in range(1,N):
                matrix[i,j] = k[j]**3 / (special.jv(q, k[j])**2 * (k[j]**2 - q**2))\
                    * (special.jv(q-1, k[j]*k[i]/S) - special.jv(q+1, k[j]*k[i]/S))
                    
    return 1/(2*np.pi * R**3) * matrix


def multiR_N(q, R, k, N, kd):
    S = kd[abs(q), N-1]
    matrix = np.zeros((N,N))
    for i in range(N):
        matrix[i,i] = k[abs(q),i]*R/S
    
    return matrix


def divR_N(q, R, k, N, kd):

    S = kd[abs(q), N-1]
    matrix = np.zeros((N,N))
    if q == 0:
        for i in range(N):
            matrix[i,i] = S/(k[abs(q),i]*R)
    
    else:
        for i in range(1,N):
            matrix[i,i] = S/(k[abs(q),i]*R)
    
        #matrix[0,0] = 1/N
    
    return matrix

def Laplacian_N(q, R, k, N):
    matrix = np.zeros((N,N))
    for i in range(N):
        matrix[i,i] = -(k[abs(q),i]/R)**2
    return matrix

def interpolation_D(r, R, n, bessel_zeros, N):
    
    jnk = bessel_zeros[abs(n), :N]
    
    W = bessel_zeros[abs(n), N]
    
    
    I = np.ndarray((len(jnk),(len(r))))
    

    for j in range(len(r)):
            
        denominator = (special.jv(n+1, jnk)*(jnk**2-(r[j]*W/R)**2))
        if abs(n) >= 0:
            denominator[0] = 1
            
        I[:,j] = 2*jnk*special.jv(n, r[j]*W/R)/denominator
    
    
    
    return I

def interpolation_N(r, R, n, k, N, kd):
    
    jnk = k[abs(n), :N]
    
    W = kd[abs(n), N-1]/R
    
    
    I = np.ndarray((len(jnk),(len(r))))
    

    for j in range(len(r)):
            
        denominator = (special.jv(n, jnk)**2 * (jnk**2-(r[j]*W)**2) * (jnk**2 - n**2))
        if abs(n) >= 0:
            denominator[0] = 1
            
        I[:,j] = 2*jnk**2 * \
            (jnk * special.jv(n+1, jnk) * special.jv(n, r[j]*W) - \
             r[j]*W * special.jv(n,jnk) * special.jv(n+1, r[j]*W))\
            /denominator
    
    
    
    return I

roots_D = np.genfromtxt('dht.csv', delimiter=',')
roots_N = np.genfromtxt('dht_neumann.csv', delimiter=',')

R = 1 #radius
N = 128 #radial modes
q = 0  #angular mode
Nr = 256 #radial sampling points

r_grid = np.linspace(0, R, Nr)


F_D =  Forward_D(q, R, roots_D, N)
B_D =  Backward_D(q, R, roots_D, N)
D_D = Derivative_D(q, R, roots_D, N)#, bessel_zeros_neu)
mR_D = multiR_D(q, R, roots_D, N)
dR_D = divR_D(q, R, roots_D, N)
dQ = dTheta_D(q, R, roots_D, N)
L_D = Laplacian_D(q, R, roots_D, N)
inter_D = interpolation_D(r_grid, R, q, roots_D, N)



F_N =  Forward_N(q, R, roots_N, N, roots_D)
B_N =  Backward_N(q, R, roots_N, N, roots_D)
L_N = Laplacian_N(q, R, roots_N, N)
D_N = Derivative_N(q, R, roots_N, N, roots_D)
mR_N = multiR_N(q, R, roots_N, N, roots_D)
dR_N = divR_N(q, R, roots_N, N, roots_D)
inter_N = interpolation_N(r_grid, R, q, roots_N, N, roots_D)

i_n = B_N @ F_N
Lap_N = B_N @ L_N @ F_N
Lap_D = B_D @ L_D @ F_D


q_test = q
k = roots_N[abs(q), 5]       #bessel root k, determine Dirichlet or Neumann

fr = special.jv(q_test, k*r_grid)  #J_q(kr) 
#fr = np.ones(len(r_grid))
interped_fr = interpolate.interp1d(r_grid, fr, kind='quadratic', fill_value='extrapolate')
#plt.figure()
#plt.plot(r_grid, fr)


root_grid_D = roots_D[abs(q), :N]*R / roots_D[abs(q), N]
root_grid_N = roots_N[abs(q), :N]*R / roots_D[abs(q), N-1]

root_grid = root_grid_D

#signal_D = interped_fr(root_grid_D) 
signal_D = special.jv(q_test, k*root_grid_D)
signal_N = special.jv(q_test, k*root_grid_N)
signal_N = np.ones(N)

result_N = B_N @ F_N @ signal_N
plt.plot(result_N)


#%%

#analytic = (special.jv(q_test-1, k*r_grid) - special.jv(q_test+1, k*r_grid)) * k/2
analytic = (signal_D*(-(k/R)**2)) #@ inter_D
output_D =   dR_D @ D_D @ F_D @ mR_D@ D_D @ F_D @ signal_D + dR_D@dR_D@dQ@ dQ @signal_D
interped_result_D  = output_D @ inter_D# + test @ inter_D)/2

#dR_N @ D_N @ F_N @ mR_N@ 
output_N =   D_N @ F_N @D_N @ F_N @ signal_N# + dR_N@dR_N@dQ@ dQ @signal_D
interped_result_N  = output_N @ inter_N

#error = np.abs(interped_result_D - analytic)
plt.figure(figsize=(10,8))
#plt.plot(error)
#plt.semilogy()
#plt.plot(root_grid, signal_D, label='Input')
plt.plot(root_grid, output_N, label='Output')
#plt.plot(root_grid, analytic, label='analytic')
#plt.plot(r_grid, interped_result_D, label='interpolated_dirich')
#plt.plot(r_grid, interped_result_N, label='interpolated_neumann')
plt.legend()








