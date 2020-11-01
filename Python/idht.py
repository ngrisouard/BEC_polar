# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 22:20:38 2020

@author: rundo
"""

'''

Inverse quasi-discrete Hankel transform of integer order n.

Input:
 H      Spectrum H(k)
 I      Integration kernel °
 K      Spectrum factors °
 R      Signal factors °

Output:
 h      Signal h(r)

°)  As computed with DHT.
'''


import numpy as np

def idht(H, I, K, R):
    
    return I*(H/K)*R