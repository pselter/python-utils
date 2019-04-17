#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 14:06:49 2018

@author: akhansen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import math
import shape_pulses as sp

#time,amp,freq,phase = sp.create_WURST(30,0.05,750,N=2)
#
#sp.plot_shape(time,amp,freq,phase)
#sp.analysize_shape(time,amp,freq,phase)


sp.calc_SHAP_par(200,10)


time,amp,freq,phase = sp.create_SHAP(40,0.05,1000,c=10,tan_k=20,create=True)
#time,amp,freq,phase = sp.create_WURST(40,0.05,1000,N=80,create=True)


sp.plot_shape(time,amp,freq,phase)
sp.analysize_shape(time,amp,freq,phase)


#fig1 = plt.figure()
#
#pulse_vec = np.zeros((len(amp),2))
#pulse_vec[:,0] = amp
#pulse_vec[:,1] = freq 
#
#print(pulse_vec)
#plt.plot(time,freq)
#plt.scatter(time,pulse_vec[:,1],c=pulse_vec[:,0],cmap='viridis')
#plt.scatter(time,pulse_vec[:,0],c=pulse_vec[:,0],cmap='seismic')
#plt.show()
#time,amp,freq,phase = sp.create_RECT(30.0,0.05,50)
#
#sp.plot_shape(time,amp,freq,phase)
#sp.analysize_shape(time,amp,freq,phase)