#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 12:55:52 2018

@author: akhansen
"""

import nmrglue as ng
import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import math

#path='/mnt/hgfs/IPC-net/sciebo/NMR-Rohdaten/DSX500/p_selt01/nmr/18-01-22-PTFE/897/pdata/998'
path='/mnt/hgfs/IPC-net/sciebo/NMR-Rohdaten/DSX500/p_selt01/nmr/18-08-15-PTFE-Rseq/898/pdata/998'
#path='/mnt/hgfs/IPC-net/sciebo/NMR-Rohdaten/DSX500/p_selt01/nmr/18-08-15-PTFE-X/899/pdata/998'
#path='C:\\Users\\philipp\\sciebo\\NMR-Rohdaten\\DSX500\\p_selt01\\nmr\\18-08-15-PTFE-Rseq\\898\\pdata\\998'


dic,data = ng.bruker.read_pdata(path)
pdata = ng.bruker.read_procs_file(path)

udic = ng.bruker.guess_udic(dic,data)
uc = ng.fileiobase.uc_from_udic(udic,0)

t_inc = 1e-6
swh = 0.5*(1/t_inc)
sw = np.linspace(-swh,swh,8000)
sw2 = np.linspace(-swh,swh,800)



x = np.arange(0,1024,1)
data_em = data*np.exp(-x/1000)
fft_data = fft(data_em.real[1:],8000)
fft_data = fftshift(fft_data.real)
fft_data2 = fft_data/fft_data.max()

binned = np.array(fft_data.real)
dat = np.mean(binned.reshape(-1, 10), axis=1)
dat = dat/dat.max()
baseline = (dat[-1]+dat[-2]+dat[-3]+dat[-4]+dat[-5]+dat[-6]+dat[-7]+dat[-8])/8
#print(dat-baseline)
dat_b=dat-baseline


#### downsample the data
start=492
end=506

fig = plt.figure(dpi=150, figsize=(10,4))
plt.plot(sw[3999:],fft_data2.real[3999:],color='b',linewidth='2',label='FFT')
plt.xlim(0,170000)

#plt.plot(sw2,dat_b+0.2,color='k',linewidth='1')
plt.scatter(sw2,dat_b,color='k',marker='s',label='downsampled')

plt.scatter(sw2[start:end],dat_b[start:end],marker='s',color='r',label='used for sim.')
plt.scatter(125782.22778473096,0.2851763073477539,marker='s',color='g',label='calibrated')
plt.xlabel(r'Nutation Frequency / kHz', fontsize=14)
plt.ylabel('Int / a.u.', fontsize=14)
test = np.linspace(0,170000,18)
plt.legend()
plt.title('B$_{1}$-Inhomogeneity', fontsize=16)
plt.xticks([0,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000,130000,140000,150000,160000,170000],['0','10','20','30','40','50','60','70','80','90','100','110','120','130','140','150','160','170'])
plt.tick_params(axis='both', which='major', labelsize=12)

#print('14 1')
for n in range((int(end)-int(start))):
    print(str(sw2[start+n]/128000)+' '+str(dat_b[start+n]))

for n in range((int(end)-int(start))):
    print(str(sw2[start+n]))

for n in range((int(end)-int(start))):
    print(str(dat_b[start+n]))

#plt.xlim(80000,170000)
#plt.show()

fig.savefig("B1_inhom.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 

#plt.close()