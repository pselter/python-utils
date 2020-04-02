# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:47:58 2020

@author: philipp
"""


import numpy as np
import TD_NMR as td
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import scipy.signal
import math



font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }


###################################
## Test of Hahn echo FID below
##
#
#test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',301)
#
#test.convdta(shifted=True,add=5)
#
#test.phase_fid(-110)
#x,y = test.return_fid()
#
#fig = plt.figure()
#plt.plot(x,y.real)
#plt.plot(x,y.imag)
#plt.show()
#

#################################
## Test of MSE FID below
##

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',220)


# x,y = test.return_fid()

fig = plt.figure(dpi=150, figsize=(5,5))
# plt.plot(x,y.imag)

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',220)

test.convdta(shifted=False,add=0)
x2,y2 = test.return_fid()

maxi = y2.imag.max()
idx = np.where(y2.imag==maxi)
print(x2[idx])
plt.plot(x2+26e-6,y2.imag)

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',221)
test.convdta(shifted=False,add=0)
x2,y2 = test.return_fid()

maxi = y2.imag.max()
idx = np.where(y2.imag==maxi)
print(x2[idx])
plt.plot(x2+20e-6,y2.imag*1.09)

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',222)
test.convdta(shifted=False,add=0)
x2,y2 = test.return_fid()

maxi = y2.imag.max()
idx = np.where(y2.imag==maxi)
print(x2[idx])
plt.plot(x2+13e-6,y2.imag*1.23)

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',223)
test.convdta(shifted=False,add=0)
x2,y2 = test.return_fid()

maxi = y2.imag.max()
idx = np.where(y2.imag==maxi)
print(x2[idx])
plt.plot(x2+6e-6,y2.imag*1.4)

test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',224)
test.convdta(shifted=False,add=0)
x2,y2 = test.return_fid()

maxi = y2.imag.max()
idx = np.where(y2.imag==maxi)
print(x2[idx])
plt.plot(x2-1e-6,y2.imag*1.6)

ax =plt.gca()
plt.xlim(-1e-6,1e-4)
# plt.ylim(-0.1,2.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)
ax.set_xlabel(r'FID / s', fontdict=font)

plt.show()
fig.savefig("NIKE-TDNMR.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# plt.close() 
