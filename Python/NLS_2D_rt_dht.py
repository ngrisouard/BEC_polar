# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 20:07:35 2020

@author: rundo
"""
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from dht import *
from ngobfft import *
import os


if os.path.isdir('./output') == False:
    os.mkdir('./output')

nti = 0             #starting iteration (loads the last file of the previous run)
ntn = 100           #final iteration of the current run


if nti == 0:        #then define all parameters

    comp_ker = 1    #=1 if the values of the kernel have to be computed

    split_NL = 0.5  #=0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
    
    if split_NL == 1:
        V = 0
    else:
        V = -0.5
    
    R = 10          #radius of the cup NORMALIZED BY THE HEALING LENGTH
    
############## Time and spatial coordinates ##############    
    nt = ntn
    dt = 5e-3
    T = nt * dt
    
    ppskip = 40     #interval between two written outputs
    
    if nti > 0:
        comp_ker = 0
        
    nth = 256       #angular points
    dth = 2 * np.pi / nth
    theta = np.arange(0, 2*np.pi, dth)
    
    nr = 256        #radial points
    dr = R / (nr - 1)
    r = np.arange(0, R+dr, dr)
    interp_sch = 'spline'
    
    jmodes = 128    #number of modes
    
    Thet, Rad = np.meshgrid(theta, r)
    
############### Initial condition #########################

    r1 = 0          #center of vortex (radius)
    theta1 = 0      #center of vortex (angle)
    circ1 = 1       #vortex circulation

    #r2 = 3*R/8    #center of vortex (radius)
    #theta2 = np.pi   #center of vortex (angle)
    #circ2 = -1    #vortex circulation

    #r3 = 3*R/8    #center of vortex (radius)
    #theta3 = 0    #center of vortex (angle)
    #circ3 = 1     #vortex circulation

    wf = np.tanh((R-Rad)/np.sqrt(2)) * np.exp(7j*Thet) \
        * (Rad * np.exp(1j*circ1*Thet) - r1 * np.exp(1j*circ1*theta1)) \
            / np.sqrt(Rad**2 + r1**2 - 2*r1*Rad*np.cos(circ1*(Thet - theta1))+1)
            
    if 'r2' in locals():
        wf = wf * (Rad*np.exp(1j*circ2*Thet) - r2*np.exp(1j*circ2*theta2)) \
            / np.sqrt(Rad**2 + r2**2 - 2*r2*Rad*np.cos(circ2*(Thet-theta2))+1)
            
    if 'r3' in locals():
        wf = wf * (Rad*np.exp(1j*Thet) - r3*np.exp(1j*theta3)) \
            / np.sqrt(Rad**2 + r3**2 - 2*r3*Rad*np.cos(circ3*(Thet-theta3))+1)
    
    
    #wf[5][5] = np.inf
    check_inf = np.zeros((256,256), int)
    np.isinf(wf, check_inf)
    
    wf = ma.masked_array(data = wf, mask = check_inf,
                    fill_value = 0, dtype = complex)
    
    Theti = np.zeros((Thet.shape[0], Thet.shape[1]+1), float)
    Theti[:Thet.shape[0], :Thet.shape[1]] = Thet
    
    Radi = np.zeros((Rad.shape[0], Rad.shape[1]+1), float)
    Radi[:Rad.shape[0], :Rad.shape[1]] = Rad
    
    wfi = np.zeros((wf.shape[0], wf.shape[1]+1), complex)
    wfi[:wf.shape[0], :wf.shape[1]] = wf
    
    for i in range(Radi.shape[0]):
        Radi[i][-1] = Rad[i][0]
      
    for i in range(wfi.shape[0]):
        wfi[i][-1] = wf[i][0]
        
    for i in range(Theti.shape[0]):
        Theti[i][-1] = Thet[i][0]
        
    X = Radi * np.cos(Theti)
    Y = Radi * np.sin(Theti)
    
    Z = np.abs(wfi)**2
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.set_aspect('equal')
    ax.plot_surface(X, Y, Z, cmap=cm.jet)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(10, -10)
    #ax.set_zlim(0,1)
    
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    plt.grid()
    plt.show()
    
        


    
    


