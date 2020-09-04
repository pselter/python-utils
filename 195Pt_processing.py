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

font1 = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }



# Add your data path here
path = 'C:/Users/pselt/OneDrive/old_sciebo_share/NMR-Rohdaten/NMR300/p_selt01/nmr/LOH-64&65'
# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'

# define the gaussian bell
# gauss = cpmg.create_gauss(npoints=402,width=20,center=170)
# gauss = cpmg.create_lorentzian(npoints=402,width=50,center=170)
# gauss = cpmg.create_rectangle(npoints=402,width=20,center=170)
gauss = cpmg.create_tukey(npoints=402,width=150,center=170,alpha=0.5)
# gauss = cpmg.create_cheby(npoints=402,width=50,center=170)
# gauss = cpmg.create_flattop(npoints=402,width=150,center=170)
# print(len(gauss))
# gauss = cpmg.create_sinc(npoints=352,width=20,center=170)
# concatenated the gaussian bells
multigauss = cpmg.create_rp_wdw(0,gauss,15,0)


test1 = cpmg.cpmg_data(path,expno=20)
# test1 = cpmg.cpmg_data(path,expno=301)


time, data = test1.return_fid()
test1.apply_wdw(multigauss)

time, data2 = test1.return_fid()
freq, fftdata = test1.return_fft(64*1024,32*1024)

test1.reload_fid()
time, data3 = test1.return_fid()
freq, fftdata2 = test1.return_fft(64*1024,32*1024)




# ###############################################################################
# fig = plt.figure(figsize=(5,4),dpi=200)
# ax=plt.gca()

# plt.plot(freq, np.abs(fftdata2)/10000,linewidth=0.25,color='k')
# plt.xlim(freq[-1],freq[0])

# plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

# ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# # ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)

# plt.show()
# # fig.savefig("CPMG_spec_1.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
# ###############################################################################


# ###############################################################################
# fig = plt.figure(figsize=(5,4),dpi=200)
# fig=plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=2)
# ax=plt.gca()
# plt.xlim(-10,750)
# plt.ylim(-0.5,10)
# plt.xlabel('FID points', fontdict=font1)

# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=120)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(multigauss*0.8+8.8,color='gray',linewidth=1)
# plt.plot(gauss*0.8+8.8,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+8.8,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5+8.2,color='navy',linewidth=0.5)



# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=140)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(multigauss*0.8+6.8,color='gray',linewidth=1)
# plt.plot(gauss*0.8+6.8,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+6.8,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5+6.2,color='navy',linewidth=0.5)


# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=160)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(multigauss*0.8+4.8,color='gray',linewidth=1)
# plt.plot(gauss*0.8+4.8,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+4.8,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5+4.2,color='navy',linewidth=0.5)


# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=180)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(multigauss*0.8+2.8,color='gray',linewidth=1)
# plt.plot(gauss*0.8+2.8,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+2.8,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5+2.2,color='navy',linewidth=0.5)



# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=200)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(multigauss*0.8+0.8,color='gray',linewidth=1)
# plt.plot(gauss*0.8+0.8,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+0.8,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5+0.2,color='navy',linewidth=0.5)


# plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)



# fig=plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)
# ax=plt.gca()


# ###############
# ## spec 1

# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=120)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+10,color='k',linewidth=0.05)



# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=140)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+9,color='red',linewidth=0.05)



# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=160)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+8,color='blue',linewidth=0.05)


# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=180)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+7,color='blue',linewidth=0.05)


# test1.reload_fid()
# gauss = cpmg.create_gauss(npoints=402,width=15,center=200)
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)
# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# plt.plot(freq,np.abs(fftdata)/np.max(np.abs(fftdata))+6,color='blue',linewidth=0.05)

# plt.xlim(freq[-1],freq[0])

# # plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
# # plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

# # ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# # # ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

# # ax.spines['top'].set_visible(False)
# # ax.spines['right'].set_visible(False)
# # ax.spines['left'].set_visible(False)

# plt.tight_layout()
# plt.show()
# # fig.savefig("processing2.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
# ###############################################################################






###############################################################################
fig = plt.figure(figsize=(5,4),dpi=200)
ax=plt.gca()
plt.xlim(-100,2000)
plt.xlabel('FID points', fontdict=font1)

plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

plt.plot(gauss*0.8+1.8,color='k',linewidth=1)
plt.plot(multigauss*0.8+1.1,color='gray',linewidth=1)
plt.plot(gauss*0.8+1.1,color='tab:red',linewidth=1)
plt.plot(data/data.max()*0.5+1.0,color='k',linewidth=0.5)
plt.plot(data2/data2.max()*0.5,color='navy',linewidth=0.5)
plt.tight_layout()

plt.show()
fig.savefig("processing_tukey.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################

###############################################################################
fig = plt.figure(figsize=(5,4),dpi=200)
ax=plt.gca()

plt.plot(freq, np.abs(fftdata)/10000,linewidth=0.25,color='navy')
plt.xlim(freq[-1],freq[0])

plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

ax.set_xlabel(r'$\nu$ '+'($^{195}$Pt) / Hz', fontdict=font1)
# ax.set_ylabel(r'Intensity / ppm', fontdict=font1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.show()
fig.savefig("CPMG_spec_tukey_single.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
###############################################################################








# fig = plt.figure()
# sumtest = np.zeros(256*1024)
# for k in range(0,200,10):
#     gauss = cpmg.create_flattop(npoints=402,width=200,center=100+k)
#     # gauss = cpmg.create_sinc(npoints=402,width=50,center=75+k)
#     # concatenated the gaussian bells
#     multigauss = cpmg.create_rp_wdw(0,gauss,150,0)
#     test1.apply_wdw(multigauss)
#     freq, fftdata = test1.return_fft(256*1024)
#     sumtest = np.fmax(np.abs(sumtest),np.abs(fftdata))
#     plt.plot(freq,np.abs(fftdata))
#     test1.reload_fid()
# # plt.plot(freq,sumtest)
# plt.xlim(freq[-1],freq[0])



# plt.show()
# # plt.close()


# # ##############################################################################
# ##############################################################################

# gauss = cpmg.create_gauss(npoints=402,width=20,center=170)
# # gauss = cpmg.create_lorentzian(npoints=402,width=5,center=170)
# # gauss = cpmg.create_rectangle(npoints=402,width=100,center=170)
# # gauss = cpmg.create_cheby(npoints=402,width=50,center=170)
# # gauss = cpmg.create_flattop(npoints=402,width=50,center=170)
# # print(len(gauss))
# # gauss = cpmg.create_sinc(npoints=352,width=20,center=170)
# # concatenated the gaussian bells
# multigauss = cpmg.create_rp_wdw(0,gauss,150,0)


# time, data = test1.return_fid()
# test1.apply_wdw(multigauss)

# time, data2 = test1.return_fid()
# freq, fftdata = test1.return_fft(64*1024,32*1024)

# test1.reload_fid()


# ###############################################################################
# fig = plt.figure(figsize=(5,4),dpi=200)
# ax=plt.gca()
# plt.xlim(-100,2000)
# plt.xlabel('FID points', fontdict=font1)

# plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

# plt.plot(gauss*0.8+1.8,color='k',linewidth=1)
# plt.plot(multigauss*0.8+1.1,color='gray',linewidth=1)
# plt.plot(gauss*0.8+1.1,color='tab:red',linewidth=1)
# plt.plot(data/data.max()*0.5+1.0,color='k',linewidth=0.5)
# plt.plot(data2/data2.max()*0.5,color='navy',linewidth=0.5)
# plt.tight_layout()

# plt.show()
# # fig.savefig("processing2.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
# ###############################################################################
