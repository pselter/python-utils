# -*- coding: utf-8 -*-
"""
Created on Sun May  3 12:14:39 2020

@author: pselt
"""


# -*- coding: utf-8 -*-
"""
Created on Sun May 19 16:35:14 2019

@author: Selter
"""

import numpy as np
from matplotlib import pyplot as plt

import sys
#sys.path.append('E:/OneDrive/projects/python-utils')
sys.path.append('C:/Users/pselt/OneDrive/projects/python-utils')
import brukerplot as bruk



font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 10,
        }


font2 = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 8,
        }






# ###############################################################################
# # 
# #  Pt-salt reference
# #


# main_fig = plt.figure(figsize=(6,4),dpi=300)


# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'
# spec = bruk.bruker1d(path, expno=20)
# ax1 = plt.gca()

# x,y = spec.plot1d_norm()

# # plt.vlines(3200.0,0,0.5, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-3700.0,0,1.1, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-5950.0,0,0.6, linestyle='--', color='silver', linewidth=1) # vertical lines

# # plt.text(2350,0.5,'$\delta_{11}$: \n3200 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-3550,1.2,'$\delta_{22}$: \n-3700 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-5700,0.65,'$\delta_{33}$: \n-5950 ppm',horizontalalignment='center', fontsize=8)

# #plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# # plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)

# plt.plot(x,y,'-',label='LOH102',color='k',linewidth=0.2)
# plt.xlim(4000,-7400)
# plt.ylim(-0.1,1.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

# plt.show()
# plt.savefig('pt_wcpmg_reference.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# ###############################################################################





# ###############################################################################
# # 
# #  Pt-salt reference
# #

# main_fig = plt.figure(figsize=(6,4),dpi=300)


# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'
# spec = bruk.bruker1d(path, expno=301)
# ax1 = plt.gca()
# ax1 = plt.gca()

# x,y = spec.plot1d_norm()

# # plt.vlines(3200.0,0,0.5, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-3700.0,0,1.1, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-5950.0,0,0.6, linestyle='--', color='silver', linewidth=1) # vertical lines

# # plt.text(2350,0.5,'$\delta_{11}$: \n3200 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-3550,1.2,'$\delta_{22}$: \n-3700 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-5700,0.65,'$\delta_{33}$: \n-5950 ppm',horizontalalignment='center', fontsize=8)

# #plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# # plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)

# plt.plot(x,y,'-',label='LOH102',color='k',linewidth=0.2)
# plt.xlim(2000,-10000)
# plt.ylim(-1.0,10.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

# plt.show()
# plt.savefig('pt_wcpmg_5wt_sample.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# ###############################################################################




# ##############################################################################
# # 
# #  Pt-salt reference
# #


# main_fig = plt.figure(figsize=(6,4),dpi=300)


# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'
# # spec = bruk.bruker1d(path, expno=20)
# ax1 = plt.gca()


# # plt.vlines(3200.0,0,0.5, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-3700.0,0,1.1, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-5950.0,0,0.6, linestyle='--', color='silver', linewidth=1) # vertical lines

# # plt.text(2350,0.5,'$\delta_{11}$: \n3200 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-3550,1.2,'$\delta_{22}$: \n-3700 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-5700,0.65,'$\delta_{33}$: \n-5950 ppm',horizontalalignment='center', fontsize=8)

# #plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# # plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)


# spec = bruk.bruker1d(path, expno=10)
# x1,y1 = spec.plot1d()
# plt.plot(x1,y1/y1.max()+2.0,'-',color='navy',linewidth=0.2)

# spec = bruk.bruker1d(path, expno=11)
# x1,y2 = spec.plot1d()
# plt.plot(x1,y2/y1.max()+1.5,'-',color='darkgreen',linewidth=0.2)

# spec = bruk.bruker1d(path, expno=12)
# x1,y3 = spec.plot1d()
# plt.plot(x1,y3/y1.max()+1.0,'-',color='firebrick',linewidth=0.2)

# spec = bruk.bruker1d(path, expno=13)
# x1,y4 = spec.plot1d()
# plt.plot(x1,y4/y1.max()+0.5,'-',color='darkorange',linewidth=0.2)

# spec = bruk.bruker1d(path, expno=14)
# x1,y5 = spec.plot1d()
# plt.plot(x1,y5/y1.max(),'-',color='k',linewidth=0.2)


# plt.xlim(4000,-8000)
# # plt.ylim(-0.1,1.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

# plt.show()
# plt.savefig('pt_wcpmg_bandwidth.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# ###############################################################################



# ###############################################################################
# # 
# #  Pt-salt reference
# #

# main_fig = plt.figure(figsize=(6,4),dpi=300)


# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'
# spec = bruk.bruker1d(path, expno=7)
# ax1 = plt.gca()
# ax1 = plt.gca()

# x,y = spec.plot1d_norm()

# # plt.vlines(3200.0,0,0.5, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-3700.0,0,1.1, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-5950.0,0,0.6, linestyle='--', color='silver', linewidth=1) # vertical lines

# # plt.text(2350,0.5,'$\delta_{11}$: \n3200 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-3550,1.2,'$\delta_{22}$: \n-3700 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-5700,0.65,'$\delta_{33}$: \n-5950 ppm',horizontalalignment='center', fontsize=8)

# #plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# # plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)

# plt.plot(x,y,'-',label='LOH102',color='k',linewidth=0.2)
# plt.xlim(3000,-7000)
# plt.ylim(-0.1,1.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

# plt.show()
# plt.savefig('pt_wcpmg_ref_750_kHz.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# ###############################################################################


# ###############################################################################
# # 
# #  Pt-salt reference
# #

# main_fig = plt.figure(figsize=(6,4),dpi=300)


# path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/400.EUSMI-PS-200109'
# spec = bruk.bruker1d(path, expno=20,procno=1)
# ax1 = plt.gca()


# x,y = spec.plot1d_norm()

# # plt.vlines(3200.0,0,0.5, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-3700.0,0,1.1, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.vlines(-5950.0,0,0.6, linestyle='--', color='silver', linewidth=1) # vertical lines

# # plt.text(2350,0.5,'$\delta_{11}$: \n3200 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-3550,1.2,'$\delta_{22}$: \n-3700 ppm',horizontalalignment='center', fontsize=8)
# # plt.text(-5700,0.65,'$\delta_{33}$: \n-5950 ppm',horizontalalignment='center', fontsize=8)

# #plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# # plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# # plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)

# plt.plot(x,y,'-',label='LOH102',color='k',linewidth=0.2)
# plt.xlim(3000,-7000)
# plt.ylim(-0.1,1.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

# plt.show()
# plt.savefig('pt_wcpmg_ref_1000_kHz.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# ###############################################################################

############################################################
# ## 19F-1H CP plot
# #  
# test = bruk.bruker2d(dpath, expno=33)
# test.def_contours2(zoom=8,mult=1.1,n_pos_levels=64,n_neg_levels=0)
# #test.def_contours_percent(zoom=5,n_pos_levels=64,n_neg_levels=64)
# #test.def_contours_percent(zoom=5,inc=2,n_pos_levels=64,n_neg_levels=64)
# fig, ax1, ax2, ax3 = test.build2d2(dpi=150,size=(6,4),limits=((-90,-155),(30,-10)),linewidth=0.5,use_cmap=True,cmap='jet',vmin=8.5,vmax=50)
# # fig, ax1, ax2, ax3 = test.build2d(dpi=150,size=(6,4),limits=((20,-5),(-90,-155)),linewidth=0.5)

# plt.sca(ax1)
# plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=12)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelsize=12)

# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{19}$F) / ppm', fontdict=font)
# ax1.set_ylabel(r'$\delta_{iso}$'+'($^{1}$H) / ppm', fontdict=font)
# #ax1.set_xlabel('$^{19$F $\delta_{iso}$ / ppm', fontname="serif", fontsize=12)
# #ax1.set_ylabel('$^{1}$H-$^{1}$H double-quantum / ppm', fontname="serif", fontsize=12)
# ax1.yaxis.set_label_position("right")
# plt.show()
# fig.savefig("direct_CP.png",dpi=300, format='png', bbox_inches='tight',  pad_inches=0.2,transparent=True) 
# #plt.close()
# ############################################################


###############################################################################
# 
#  Pt-salt reference indirect detection
#




path = 'D:/Box Sync/Research/Lille-NMR/IR_RMN/800.EUSMI-PS-200113'
spec = bruk.bruker2d(path, expno=17,procno=1,normalize=True)


spec.def_contours2(zoom=5,mult=1.5,n_pos_levels=32,n_neg_levels=32)
fig, ax1, ax2, ax3 = spec.build2d2(dpi=150,size=(8,6),limits=((15,-5),(4000,-9000)),linewidth=0.5,use_cmap=True,cmap='jet',vmin=0.5,vmax=50)


plt.sca(ax1)
plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=8)
plt.tick_params(axis='y', which='both', labelleft=False, labelright=True, left=False,right=True, labelsize=8)

ax1.yaxis.set_label_position("right")

ax1.set_xlabel('$^{1}$H $\delta_{iso}$ / ppm', fontname="serif", fontsize=12)
ax1.set_ylabel('$^{195}$Pt / ppm', fontname="serif", fontsize=12)


#plt.vlines(-2100.0,0,20000, linestyle='--', color='silver', linewidth=1) # vertical lines
# plt.text(-1500,0.6,'$\delta_{iso}$: -2150 ppm \n$\delta_{aniso}$: 5350 ppm\n$\eta$: 0.42',horizontalalignment='center', fontsize=8,color='tab:red')

# plt.text(500,0.5,'Pt$^{II}$[NH$_{3}$]$_{4}$Cl$_{2}$',horizontalalignment='center', fontdict=font2)
x,y = spec.slice_0(1024)
plt.plot(y/20+7.5,x,'-',color='maroon',linewidth=0.5)
x,y = spec.slice_0(900)
plt.plot(y/20,x,'-',color='k',linewidth=0.5)
# plt.xlim(3000,-7000)
# plt.ylim(-0.1,1.0)


# plt.tick_params(axis='x', which='both', labelbottom=True, bottom=True, labelsize=8)
# plt.tick_params(axis='y', which='both', labelleft=False, left=False)
# plt.minorticks_on()
# ax1.set_xlabel(r'$\delta_{iso}$'+'($^{195}$Pt) / ppm', fontdict=font)
# #plt.legend()


# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['left'].set_visible(False)

plt.show()
plt.savefig('pt_indirect.png',dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
###############################################################################
