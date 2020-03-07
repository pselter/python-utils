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

# Add your data path here
path = 'C:/Users/Selter/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH-64&65'

# define the gaussian bell
# gauss = cpmg.create_gauss(npoints=402,width=2,center=150)
# gauss = cpmg.create_lorentzian(npoints=402,width=5,center=170)
# gauss = cpmg.create_rectangle(npoints=402,width=100,center=170)
gauss = cpmg.create_cheby(npoints=402,width=160,center=170)
print(len(gauss))
# gauss = cpmg.create_sinc(npoints=402,width=50,center=170)
# concatenated the gaussian bells
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
# plt.xlim(0,300)
plt.show()


fig = plt.figure()
plt.plot(freq, np.abs(fftdata))
plt.xlim(freq[-1],freq[0])
plt.show()


fig = plt.figure()
sumtest = np.zeros(64*1024)
for k in range(0,200,2):
    gauss = cpmg.create_cheby(npoints=402,width=150,center=75+k)
    # gauss = cpmg.create_sinc(npoints=402,width=10,center=75+k)
    # concatenated the gaussian bells
    multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
    test1.apply_wdw(multigauss)
    freq, fftdata = test1.return_fft(256*1024)
    # sumtest = sumtest + np.abs(fftdata)
    test1.reload_fid()
    plt.plot(freq,np.abs( fftdata))
plt.xlim(freq[-1],freq[0])
plt.show()
