# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:57:08 2020

@author: pselt
"""
import numpy as np
import cpmg_utils as cpmg
import matplotlib.pyplot as plt


font1 = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }


###############################################################################
# Path to example data
#
path = 'C:/Users/pselt/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH93'
# path = 'C:/Users/pselt/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH-64&65'
# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'

# define the gaussian bell
# gauss = cpmg.create_gauss(npoints=402,width=20,center=170)
# gauss = cpmg.create_lorentzian(npoints=402,width=5,center=170)
# gauss = cpmg.create_rectangle(npoints=402,width=100,center=170)
# gauss = cpmg.create_cheby(npoints=402,width=50,center=170)
# gauss = cpmg.create_flattop(npoints=352,width=50,center=170)
# print(len(gauss))

gauss = cpmg.create_rectangle(npoints=452,width=452,center=226)
# gauss = cpmg.create_tukey(npoints=452,width=150,center=226,alpha=0.5)

# gauss = cpmg.create_sinc(npoints=402,width=100,center=170)
# concatenated the gaussian bells
multigauss = cpmg.create_rp_wdw(0,gauss,10,0)


test1 = cpmg.cpmg_data(path,expno=10)
# test1 = cpmg.cpmg_data(path,expno=301)


time, data = test1.return_fid()
test1.apply_wdw(multigauss)

time, data2 = test1.return_fid()
freq, fftdata5 = test1.return_fft(256*1024,128*1024)

test1.reload_fid()
time, data3 = test1.return_fid()
freq, fftdata2 = test1.return_fft(256*1024,128*1024)



###############################################################################
#
# Original Spectrum
#
###################
fig = plt.figure(figsize=(5,3),dpi=200)
ax=plt.gca()

plt.plot(freq, np.abs(fftdata5)/10000,linewidth=0.1,color='k')
plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.show()
# fig.savefig("CPMG_spec_org_TDEFF_100echoes.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###################
###############################################################################



###############################################################################
#
# New processing
#
###################
fig = plt.figure(figsize=(5,3),dpi=200)
ax=plt.gca()
sumtest = np.zeros(256*1024)


for k in range(0,200,1):
    gauss = cpmg.create_sinc(npoints=452,width=15,center=126+k)
    # gauss = cpmg.create_flattop(npoints=452,width=100,center=126+k)
    # gauss = cpmg.create_rectangle(npoints=452,width=100,center=126+k)
    # gauss = cpmg.create_tukey(npoints=352,width=100,center=126+k,alpha=0.5)
    # gauss = cpmg.create_cheby(npoints=452,width=75,center=126+k)
    # gauss = cpmg.create_lorentzian(npoints=452,width=50,center=126+k)
    # gauss = cpmg.create_gauss(npoints=352,width=50,center=76+k)
    # concatenated the gaussian bells
    multigauss = cpmg.create_rp_wdw(0,gauss,10,0)
    test1.apply_wdw(multigauss)
    freq, fftdata = test1.return_fft(256*1024)
    # sumtest = sumtest+np.abs(fftdata)
    sumtest = np.fmax(np.abs(sumtest),np.abs(fftdata))
    # plt.plot(freq,np.abs(fftdata))
    test1.reload_fid()
plt.plot(freq, np.abs(fftdata5)/1000000,linewidth=0.2,color='r')
plt.plot(freq,sumtest/1000000,linewidth=0.2,color='k')


plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.show()
# fig.savefig("CPMG_spec_sinc15_td100.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################


