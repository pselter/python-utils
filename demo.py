# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:48:58 2020

@author: pselt
"""


import numpy as np
import nmrglue as ng
import cpmg_utils as cpmg
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift

font1 = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }



# Add your data path here
path = 'C:/Users/pselt/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH-64&65'
# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'


test1 = cpmg.cpmg_data(path,expno=20)
width = 15




###############################################################################
fig = plt.figure(figsize=(5,4),dpi=200)
plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=2)
ax=plt.gca()
plt.xlim(-10,750)
plt.ylim(-0.5,10)
plt.xlabel('FID points', fontdict=font1)

test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+8.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+8.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+8.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+8.2,color='navy',linewidth=0.5)



test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+6.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+6.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+6.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+6.2,color='navy',linewidth=0.5)


test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+4.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+4.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+4.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+4.2,color='navy',linewidth=0.5)


test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+2.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+2.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+2.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+2.2,color='navy',linewidth=0.5)



test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+0.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+0.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+0.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+0.2,color='navy',linewidth=0.5)


plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)



plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
ax=plt.gca()


###############
## spec 1

test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+10,color='k',linewidth=0.05)



test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+9,color='red',linewidth=0.05)



test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+8,color='blue',linewidth=0.05)


test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+7,color='blue',linewidth=0.05)


test1.reload_fid()
gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+6,color='blue',linewidth=0.05)

plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
plt.show()
fig.savefig("processing_gauss.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################










width = 20





###############################################################################
fig = plt.figure(figsize=(5,4),dpi=200)
plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=2)
ax=plt.gca()
plt.xlim(-10,750)
plt.ylim(-0.5,10)
plt.xlabel('FID points', fontdict=font1)

test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
gauss = cpmg.create_sinc(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+8.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+8.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+8.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+8.2,color='navy',linewidth=0.5)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
gauss = cpmg.create_sinc(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+6.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+6.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+6.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+6.2,color='navy',linewidth=0.5)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
gauss = cpmg.create_sinc(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+4.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+4.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+4.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+4.2,color='navy',linewidth=0.5)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
gauss = cpmg.create_sinc(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+2.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+2.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+2.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+2.2,color='navy',linewidth=0.5)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
gauss = cpmg.create_sinc(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+0.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+0.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+0.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+0.2,color='navy',linewidth=0.5)


plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)



plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
ax=plt.gca()


###############
## spec 1

test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
gauss = cpmg.create_sinc(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+10,color='k',linewidth=0.05)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
gauss = cpmg.create_sinc(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+9,color='red',linewidth=0.05)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
gauss = cpmg.create_sinc(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+8,color='blue',linewidth=0.05)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
gauss = cpmg.create_sinc(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+7,color='blue',linewidth=0.05)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
gauss = cpmg.create_sinc(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+6,color='blue',linewidth=0.05)

plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
plt.show()
fig.savefig("processing_sinc.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################







width = 20





###############################################################################
fig = plt.figure(figsize=(5,4),dpi=200)
plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=2)
ax=plt.gca()
plt.xlim(-10,750)
plt.ylim(-0.5,10)
plt.xlabel('FID points', fontdict=font1)

test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=120)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+8.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+8.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+8.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+8.2,color='navy',linewidth=0.5)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=140)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+6.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+6.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+6.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+6.2,color='navy',linewidth=0.5)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=160)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+4.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+4.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+4.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+4.2,color='navy',linewidth=0.5)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=180)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+2.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+2.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+2.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+2.2,color='navy',linewidth=0.5)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=200)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(multigauss*0.8+0.8,color='gray',linewidth=1)
plt.plot(gauss*0.8+0.8,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+0.8,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5+0.2,color='navy',linewidth=0.5)


plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)



plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
ax=plt.gca()


###############
## spec 1

test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=120)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=120)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=120)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+10,color='k',linewidth=0.05)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=140)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=140)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=140)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+9,color='red',linewidth=0.05)



test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=160)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=160)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=160)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+8,color='blue',linewidth=0.05)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=180)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=180)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=180)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+7,color='blue',linewidth=0.05)


test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=width,center=200)
gauss = cpmg.create_rectangle(npoints=402,width=width,center=200)
# gauss = cpmg.create_sinc(npoints=402,width=width,center=200)
multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
time, data = test1.return_fid()
test1.apply_wdw(multigauss)
time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+6,color='blue',linewidth=0.05)

plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.tight_layout()
plt.show()
fig.savefig("processing_rectangle.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################
