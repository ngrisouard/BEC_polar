# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 13:09:32 2021

@author: rundo
"""
import numpy as np


def compute_energy(wf, dr, dtheta, Rad):
    
    
    #wf = np.abs(wf)
    dwfdtheta = np.zeros(np.shape(wf))
    dwfdr = np.zeros(np.shape(wf))
    
    dwfdr, dwfdtheta = np.gradient(wf)
    
    for i in range(np.shape(dwfdtheta)[1]):
        
        dwfdtheta[i,:] = (np.roll(wf[i,:], -1) - wf[i,:])/dtheta
        
    for i in range(np.shape(dwfdr)[0]):
        
        dwfdr[:,i] = (np.roll(wf[:,i], -1) - wf[:,i])/dr
        
    dwfdr[0] = 0
    dwfdr[-1] = 0
    
    
    temp1 = np.zeros(np.shape(wf))
    temp2 = np.zeros(np.shape(wf))
    temp3 = np.zeros(np.shape(wf))
    temp4 = np.zeros(np.shape(wf))
    
    temp1[1:,:] = Rad[1:]*np.abs(dwfdr[1:,:])**2
    temp2[1:,:] = (np.abs(dwfdtheta[1:,:])**2)/Rad[1:]
    temp3[1:,:] = Rad[1:]*np.abs(wf[1:,:])**2
    temp4[1:,:] = 0.5*Rad[1:]*np.abs(wf[1:,:])**4
    
    I = 0.5*dr*dtheta*np.sum(temp1+temp2-temp3+temp4)

    return I    

def compute_L(wf, dr, dtheta, Rad):
    dwfdtheta = np.zeros(np.shape(wf))
    dwfdr = np.zeros(np.shape(wf))
    
    dwfdr, dwfdtheta = np.gradient(wf)
    
    for i in range(np.shape(dwfdtheta)[1]):
        
        dwfdtheta[i,:] = (np.roll(wf[i,:], -1) - wf[i,:])/dtheta
        
    for i in range(np.shape(dwfdr)[0]):
        
        dwfdr[:,i] = (np.roll(wf[:,i], -1) - wf[:,i])/dr
        
    dwfdr[0] = 0
    dwfdr[-1] = 0
    
    L = np.sum(dwfdtheta*np.conj(wf)*Rad)*dr*dtheta
    
    return L
    
    
              
        
    