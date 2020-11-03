# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:02:58 2020

@author: rundo
"""

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import io
from dht import *
from ngobfft import *
import os


misc = io.loadmat('output/misc.mat')
wf = io.loadmat('output/psi_00000000.mat')