
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
from scipy import io
from scipy import interpolate
from dht_disk import *
from dht import *
from ngobfft import *
from idht import *
import os


if os.path.isdir('./output') == False:
    os.mkdir('./output')
    os.mkdir('./output/kernal')
    


nti = 0             #starting iteration (loads the last file of the previous run)
ntn = 100           #final iteration of the current run


if nti == 0:        #then define all parameters

    comp_ker = 1    #=1 if the values of the kernel have to be computed

    split_NL = 0.5  #=0.5 for NLS, =1 for LS, =1/3 for NLS with dissipation
    
    if split_NL == 1:
        V = 0
    else:
        V = -0.5
    
    R = 1.25          #radius of the cup NORMALIZED BY THE HEALING LENGTH
    R0 = 1
    
############## Time and spatial coordinates ##############    
    nt = ntn
    dt = 2e-1
    T = nt * dt
    
    ppskip = 40     #interval between two written outputs
    
    if nti > 0:
        comp_ker = 0
        
    nth = 256       #angular points
    dth = 2 * np.pi / nth
    theta = np.arange(0, 2*np.pi, dth)
    
    nr = 256        #radial points
    
    dr = (R-R0) / (nr - 1)
    #dr = R / (nr - 1)
    
    r = np.arange(R0, R+dr, dr)
    #r = np.arange(0, R+dr, dr)
    
    interp_sch = 'quadratic'
    
    jmodes = 128    #number of modes
    
    Thet, Rad = np.meshgrid(theta, r)
    
############### Initial condition #########################

    r1 = 3*R/4         #center of vortex (radius)
    theta1 = 0      #center of vortex (angle)
    circ1 = 1       #vortex circulation

    #r2 = 0    #center of vortex (radius)
    #theta2 = 0   #center of vortex (angle)
    #circ2 = 0    #vortex circulation

    #r3 = 3*R/8    #center of vortex (radius)
    #theta3 = 0    #center of vortex (angle)
    #circ3 = 1     #vortex circulation
    
    #* np.tanh((Rad-R0)/np.sqrt(2))
    
    wf = np.tanh((R-Rad)/np.sqrt(2)) * np.exp(0j*Thet) * np.tanh((Rad-R0)/np.sqrt(2))\
    #    * (Rad * np.exp(1j*circ1*Thet) - r1 * np.exp(1j*circ1*theta1)) \
    #        / np.sqrt(Rad**2 + r1**2 - 2*r1*Rad*np.cos(circ1*(Thet - theta1))+1)
            
    if 'r2' in locals():
        wf = wf * (Rad*np.exp(1j*circ2*Thet) - r2*np.exp(1j*circ2*theta2)) \
            / np.sqrt(Rad**2 + r2**2 - 2*r2*Rad*np.cos(circ2*(Thet-theta2))+1)
            
    if 'r3' in locals():
        wf = wf * (Rad*np.exp(1j*Thet) - r3*np.exp(1j*theta3)) \
            / np.sqrt(Rad**2 + r3**2 - 2*r3*Rad*np.cos(circ3*(Thet-theta3))+1)
    '''
    for i in range(nr):
        if i <= int(nr*R0//R):
            wf[i, :] = 0
    '''
    
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
    #Z = np.angle(wfi)
    z_max = np.max(Z)
    z_min = np.min(Z)


    plt.figure(figsize = (10, 8))
    plt.pcolor(X, Y, Z, cmap=cm.jet)#, vmax = z_max, vmin = z_min)
    plt.title('t=0')
    plt.colorbar()
    plt.show()
    
    ##################### Compute the integration kernal ###########
    comp_ker = 1
    if comp_ker == 1:

        bessel_zeros = np.genfromtxt('dht_ring.csv', delimiter=',')

        for ii in range(-nth//2+1, nth//2+1):

            #H,kk,rr,I,KK,RR, h = dht_disk([],R0,R,bessel_zeros,jmodes,ii)
            H,kk,rr,I,KK,RR, h = dht_disk([],R0,R,bessel_zeros,jmodes,ii)

            save_dict = {'kk':kk, 'rr':rr, 'KK':KK, 'RR':RR, 'I':I}
            io.savemat('./output/kernal/ker_%i.mat'%ii, save_dict)
            
        
    t = 0
    nameout = 'output/psi_00000000.mat'
    io.savemat(nameout, {'wf':wf, 't':t})
    npi = 8          # # of digits following psi_
    
    misc = {'r':r, 'theta':theta, 'T':T, 'nt':nt, 'dt':dt, 'nth':nth, 
                 'dth':dth, 'nr':nr, 'dr':dr, 'R':R, 'ppskip':ppskip, 'Rad':Rad,
                 'Thet':Thet, 'npi':npi, 'V':V, 'comp_ker':comp_ker, 'interp_sch':
                     interp_sch, 'wf':wf, 't':t, 'split_NL':split_NL, 
                     'r1':r1, 'theta1':theta1, 'circ1':circ1}
    if 'r3' in locals():
        misc.update({'r2':r2, 'theta2':theta2, 'circ2':circ2, 'theta3':theta3,
                          'circ3':circ3})
    elif 'r2' in locals():
        misc.update({'r2':r2, 'theta2':theta2, 'circ2':circ2})
    
    io.savemat('output/misc.mat', misc)

    
############### then reload all the parameters and last iteration
else:
    misc = io.loadmat('output/misc.mat')
    nt = misc['ntn']
    T = misc['dt'] * misc['nt']
    misc.update({'T':T, 'nt':nt})
    io.savemat('output/misc.mat', misc)
    
    psi = io.loadmat(nameout)
    wf = psi['wf']
    t = psi['t'][0]
    
#%% Temporal loop
nti = 0
nr = 256
nth = 256
pp = nti + 1
bessel_zeros = 0 #np.genfromtxt('dht.csv', delimiter=',')
deriv_thet = 1j*(np.outer(np.ones(nr), np.arange(-nth//2+1, nth//2+1)))


 



for t in np.arange((nti+1)*dt, T+dt, dt):
    # 1st step: the linear step, computation of the laplacian part.
    
    kt = k_of_x(Thet[0]) 
    fr = obfft(Thet[0], wf, -1) #fft to isolate every angular mode
    
    if abs((pp-1)/ppskip - np.floor((pp-1)/ppskip)) <= 1e-14:
        r2 = x_of_k(kt)    
        dwfdt = obifft(kt, deriv_thet*fr, -1)




    for ii in range(-nth//2 +1, nth//2 +1):
    
        #ii = -11
        ker = io.loadmat('output/kernal/ker_%i.mat'%(ii))
        I = ker['I']
        kk = ker['kk']
        KK = ker['KK']
        rr = ker['rr']
        RR = ker['RR']
        
        #if abs(ii) > 0:
        #    I[0, 0] = 0
    
        #forward dht
        Fr_func = interpolate.interp1d(r, fr[:, ii-1+nth//2], kind=interp_sch, fill_value='extrapolate')
        Fr = Fr_func(rr)[0]

        adht,_,_,_,_,_,_ = dht(Fr, RR, bessel_zeros, KK, I)

        
        #inverse dht
        
        
        
        Fr = idht(adht*np.exp(-1j*0.5*(kk**2)*dt), I, KK, RR)
        #Fr = idht(adht, I, KK, RR)
        fr_func = interpolate.interp1d(rr[0], Fr[0], kind=interp_sch, fill_value='extrapolate')
        fr[:, ii-1+nth//2] = fr_func(r)
    
    r2 = x_of_k(kt)
    wf = obifft(kt, fr,-1)

    if abs((pp-1)/ppskip - np.floor((pp-1)/ppskip)) <= 1e-14:
    

        psi = io.loadmat(nameout)
        psi.update({'dwfdt':dwfdt})
        #io.savemat(nameout, psi)
    

    #2nd step NL part
    if split_NL == 0.5 or split_NL == 0 or abs(split_NL-1/3) < 1e-13:
        wf = wf * np.exp(-1j*(V+0.5*np.abs(wf)**2)*dt)
    
    #3rd step dssipation not yet implemented
    ###
    '''
    for i in range(nr):
        if i <= 128:
            wf[i, :] = 0
    '''
    
    wfi = np.zeros((wf.shape[0], wf.shape[1]+1), complex)
    wfi[:wf.shape[0], :wf.shape[1]] = wf

    for i in range(wfi.shape[0]):
        wfi[i][-1] = wf[i][0]

    #Z = np.angle(wfi)
    Z = np.abs(wfi)**2
    plt.figure(figsize = (10, 8))
    plt.pcolor(X, Y, Z, cmap=cm.jet)
    plt.title('t=%f'%t)
    plt.colorbar()
    plt.pause(0.05)
    plt.show()
    
    pp += 1
