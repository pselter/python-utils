# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import brukerplot as bruk
import nmrglue as ng
from scipy.optimize import curve_fit



###############################################################################
# Define fit_function that can be used for fitting
# most common are single_expo, dual_expo, and triple_expo
#

#-------------------------
# define fitfunction
def spec_dense(x,k,tc):
    f=(1+((x**2)*(tc**2)))
    return k*tc/f
#------------------------

#-------------------------
# define fitfunction
def l_peak(x,x0,w,I):
    f=(x0-x)/(w*0.5)
    return I/(1+f**2)
#------------------------
    
#-------------------------
# define fitfunction
def g_peak(x,x0,w,I):
    g=(x0-x)/(w*0.5)
    return I*np.exp(-1*(np.log(2)*g**2))
#------------------------

#-------------------------
# define the fitfunction
def triple_expo(x,T1r_1,I1,T1r_2,I2,T1r_3,I3):
    return I1*(np.exp(-1*x/T1r_1))+I2*(np.exp(-1*x/T1r_2))+I3*(np.exp(-1*x/T1r_3))
#-------------------------
 
#-------------------------
# define the fitfunction
def dual_expo(x,T1r_1,I1,T1r_2,I2):
    return I1*(np.exp(-1*x/T1r_1))+I2*(np.exp(-1*x/T1r_2))
#-------------------------
    
#-------------------------
# define the fitfunction
def single_expo(x,T1r_1,I1):
    return I1*(np.exp(-1*x/T1r_1))
#-------------------------

#-------------------------
# define two peak model
def gl_peak(x,x0_1,w1,I1,x0_2,w2,I2):
    return g_peak(x,x0_1,w1,I1)+l_peak(x,x0_2,w2,I2)
#------------------------



###############################################################################
    





def T1r_analysis(dpath,expno,region,procno,npoints,fit_type='single_expo'):
    vplist = np.loadtxt(dpath+'/'+str(expno)+'/vplist')
    integral = np.zeros(npoints)
    time = np.arange(0,vplist[-1],vplist[0])
    data = np.zeros((npoints,2))
    res = np.zeros((npoints,2))
    fit =  np.zeros((len(time),2))
    # integrate all spectra in the specified region
    for n in range(npoints):
        spectrum = bruk.bruker1d(dpath,expno,procno=procno+n)
        x,y = spectrum.plot1d()
        xaxis = np.array(x)
        idx0 = (np.abs(xaxis - region[0])).argmin()
        idx1 = (np.abs(xaxis - region[1])).argmin()
        for l in range(idx1-idx0):
            integral[n] += y[idx0+l]
        
        # normalize to 0-100
        # this makes fit estimation more stable
        norm_int = integral/integral.max()*100
    data[:,0] = vplist
    data[:,1] = norm_int


    # do the fitting using the specified fitting type
    
    if fit_type == 'single_expo':
        popt, pcov = curve_fit(single_expo,vplist,norm_int,maxfev=5000,p0=(1.0e-2,100))
        fitted = single_expo(time,*popt)
        residual = norm_int-single_expo(vplist,*popt)

    if fit_type == 'dual_expo':
        popt, pcov = curve_fit(dual_expo,vplist,norm_int,maxfev=5000,p0=(1.0e-2,50,1.0e-3,50))
        fitted = dual_expo(time,*popt)
        residual = norm_int-dual_expo(vplist,*popt)    
        
    fit[:,0] = time
    fit[:,1] = fitted
    res[:,0] = vplist
    res[:,1] = residual
        
    return data, fit, res, popt, pcov 