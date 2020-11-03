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

#c = genfromtxt('dht.csv', delimiter=',')

R = 10
jmodes = 128
n = 215
ii = n
N=256
h=[]





#H, kk, rr, I, KK, RR, h = dht([],R,jmodes,ii)

#save_dict = {'kk':kk, 'rr':rr, 'KK':KK, 'RR':RR, 'I':I}

#io.savemat('./output/ker.mat', save_dict)

read = io.loadmat('./output/ker_3.mat')

I = read['I']


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

    
