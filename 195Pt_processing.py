# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:58:54 2020

@author: Selter
"""

import numpy as np
import nmrglue as ng
import cpmg_utils as cpmg
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift

# path = 'C:/Users/Selter/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH102-2'
path = 'C:/Users/philipp/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH-64&65'

gauss = cpmg.create_gauss(npoints=402,width=10,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)


test1 = cpmg.cpmg_data(path,expno=20)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)




fig = plt.figure()
plt.plot(data/data.max()*0.5)
plt.plot(multigauss)
plt.plot(data2/data2.max()*0.5)
plt.show()



fig = plt.figure()
plt.plot(freq, np.abs(fftdata))
plt.xlim(freq[-1],freq[0])
plt.show()



# fig = plt.figure()
# # (0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)

# for this in (np.arange(0,150,10)):
#     gauss = cpmg.create_gauss(npoints=402,width=15,center=150+this)
#     print(100+this)
#     multigauss = cpmg.create_rp_wdw(0,gauss,100,0)
#     test1.apply_wdw(multigauss)
#     freq, data = test1.return_fft(128*1024)
#     plt.plot(freq, np.abs(data))
#     test1.reload_fid()

# plt.xlim(freq[-1],freq[0])
# plt.show()