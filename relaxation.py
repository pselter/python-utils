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
# define the fitfunction
def satrec1(x,T1_1,I1):
    return I1*(1-np.exp(-1*x/T1_1))
#-------------------------

#-------------------------
# define the fitfunction
def satrec2(x,T1_1,I1,T1_2,I2):
    return I1*(1-np.exp(-1*x/T1_1))+I2*(1-np.exp(-1*x/T1_2))
#-------------------------

#-------------------------
# define the fitfunction
def satrec3(x,T1_1,I1,T1_2,I2,T1_3,I3):
    return I1*(1-np.exp(-1*x/T1_1))+I2*(1-np.exp(-1*x/T1_2))+I3*(1-np.exp(-1*x/T1_3))
#-------------------------
 

#-------------------------
# define two peak model
def gl_peak(x,x0_1,w1,I1,x0_2,w2,I2):
    return g_peak(x,x0_1,w1,I1)+l_peak(x,x0_2,w2,I2)
#------------------------

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

###############################################################################
    

class relaxation(object):
    """General Object to handle relaxation measurements by pseudo2D
    
    expects the data being pre-processed via rowext in Topsin 
    or euqivalent treatment prior to use
    
    vdlists or vplists should be in seconds, no units just numbers
    """
    
    
    #----------------------------------
    def __init__(self, path, expno):
        """"Initialize the class
        
        nothing fancy, yet
        """
        self.path = path
        self.expno = expno
    #----------------------------------

    
    def T1rho(self,start_procno,npoints,region=(20,-10), normalize=True):
        """Function to do the integration, returns the integral and time in an array         
        written with T1_rho measurements by spin-locking in mind
        
        assumes a vplist being present
        - start_procno :the first procno with slices
        - npoints      :the number of procnos to evaluate
        - region       :the region to integrate in ppm
        - normalize    :normalizes the data to 100
        """
        
        self.vplist = np.loadtxt(self.path+'/'+str(self.expno)+'/vplist')
        self.integral = np.zeros(npoints)
        #self.time = np.arange(0,self.vplist[-1],self.vplist[0])
        self.data = np.zeros((npoints,2))
        #self.res = np.zeros((npoints,2))
        #self.fit =  np.zeros((len(self.time),2))
        
        #----------------------------------
        #Do the integration over all procnos
        for n in range(npoints):
            spectrum = bruk.bruker1d(self.path,self.expno,procno=start_procno+n)
            self.x, self.y = spectrum.plot1d()
            self.xaxis = np.array(self.x)
            self.idx0 = (np.abs(self.xaxis - region[0])).argmin()
            self.idx1 = (np.abs(self.xaxis - region[1])).argmin()
            for l in range(self.idx1-self.idx0):
                self.integral[n] += self.y[self.idx0+l]
        #----------------------------------
        
        #----------------------------------
        # Check for normalization
        if normalize == True:
            self.norm_int = self.integral/self.integral.max()*100
            print(self.vplist)
            self.data[:,0] = self.vplist
            self.data[:,1] = self.norm_int
        else:
            self.data[:,0] = self.vplist
            self.data[:,1] = self.integral
        #----------------------------------
        return self.data


    def T1_satrec(self,start_procno,npoints,region=(20,-10), normalize=True):
        """Function to do the integration, returns the integral and time in an array         
        written with T1 measurements by saturation recovery in mind
        
        assumes a vplist being present
        - start_procno :the first procno with slices
        - npoints      :the number of procnos to evaluate
        - region       :the region to integrate in ppm
        - normalize    :normalizes the data to 100
        """
        
        
        self.vdlist = np.loadtxt(self.path+'/'+str(self.expno)+'/vdlist')
        self.integral = np.zeros(npoints)
        #self.time = np.arange(0,self.vplist[-1],self.vplist[0])
        self.data = np.zeros((npoints,2))
        #self.res = np.zeros((npoints,2))
        #self.fit =  np.zeros((len(self.time),2))
        
        #----------------------------------
        #Do the integration over all procnos
        for n in range(npoints):
            spectrum = bruk.bruker1d(self.path,self.expno,procno=start_procno+n)
            self.x, self.y = spectrum.plot1d()
            self.xaxis = np.array(self.x)
            self.idx0 = (np.abs(self.xaxis - region[0])).argmin()
            self.idx1 = (np.abs(self.xaxis - region[1])).argmin()
            for l in range(self.idx1-self.idx0):
                self.integral[n] += self.y[self.idx0+l]
        #----------------------------------
        
        #----------------------------------
        # Check for normalization
        if normalize == True:
            self.norm_int = self.integral/self.integral.max()*100
            self.data[:,0] = self.vdlist
            self.data[:,1] = self.norm_int
        else:
            self.data[:,0] = self.vdlist
            self.data[:,1] = self.integral
        #----------------------------------
        return self.data         

    
    def exp_decay(self,data,n_expo=1):
        """ Function to fit an exponential decay of your data
        
        - data : array with the times and intensities
        - n_expo : 1-3, how many exponential decays to fit
        
        returns the fit, the residual, and the results
        """
        
        self.data = data
        self.npoints = len(data[:,0])
        self.time = np.arange(0,self.data[-1,0],self.data[0,0])
       
        self.res = np.zeros((self.npoints,2))
        self.fit =  np.zeros((len(self.time),2))
        
        self.errors = np.zeros(n_expo)
        self.T2 = np.zeros(n_expo)
        self.Int = np.zeros(n_expo)
        
        self.fit[:,0] = self.time
        self.res[:,0] = self.data[:,0]
        
        if n_expo == 1:
            self.popt, self.pcov = curve_fit(single_expo,self.data[:,0],self.data[:,1],maxfev=5000,p0=(1.0e-2,100))
            self.fit[:,1] = single_expo(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-single_expo(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2 = self.popt[0]
            self.Int = self.popt[1]
    
        if n_expo == 2:
            self.popt, self.pcov = curve_fit(dual_expo,self.data[:,0],self.data[:,1],maxfev=5000,p0=(1.0e-2,50,1.0e-3,50))
            self.fit[:,1] = dual_expo(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-dual_expo(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
            self.Int[1] = self.popt[3]
        
        if n_expo == 3:
            self.popt, self.pcov = curve_fit(triple_expo,self.data[:,0],self.data[:,1],maxfev=5000,p0=(5.0e-2,33,1.0e-2,33,5.0e-3,33))
            self.fit[:,1] = triple_expo(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-triple_expo(self.data[:,0],*self.popt)

            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
            self.Int[1] = self.popt[3]
            self.T2[2] = self.popt[4]
            self.Int[2] = self.popt[5]

        return self.fit, self.res, self.T2, self.Int, self.errors

    def sat_rec(self,data,n_expo=1):
        """ Function to fit saturation recovery experiments
        
        - data : array with the times and intensities
        - n_expo : 1-3, how many exponentials to use
        
        returns the fit, the residual, and the results
        """
        
        self.npoints = len(data[:,0])
        self.time = np.arange(0,self.data[-1,0],self.data[0,0])
        
        self.data = data
        self.res = np.zeros((self.npoints,2))
        self.fit =  np.zeros((len(self.time),2))
        
        self.fit[:,0] = self.time
        self.res[:,0] = self.data[:,0]
        
        self.errors = np.zeros(n_expo)
        self.T2 = np.zeros(n_expo)
        self.Int = np.zeros(n_expo)
        
        if n_expo == 1:
            self.popt, self.pcov = curve_fit(satrec1,self.data[:,0],self.data[:,1],maxfev=5000,p0=(1.0e-2,100))
            self.fit[:,1] = satrec1(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-satrec1(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2 = self.popt[0]
            self.Int = self.popt[1]
    
        if n_expo == 2:
            self.popt, self.pcov = curve_fit(satrec2,self.data[:,0],self.data[:,1],maxfev=5000,p0=(1.0e-2,50,1.0e-3,50))
            self.fit[:,1] = satrec2(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-satrec2(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
            self.Int[1] = self.popt[3]
        
        if n_expo == 3:
            self.popt, self.pcov = curve_fit(satrec3,self.data[:,0],self.data[:,1],maxfev=5000,p0=(5.0e-2,33,1.0e-2,33,5.0e-3,33))
            self.fit[:,1] = satrec3(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-satrec3(self.data[:,0],*self.popt)
    
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
            self.Int[1] = self.popt[3]
            self.T2[2] = self.popt[4]
            self.Int[2] = self.popt[5]
    
        return self.fit, self.res, self.T2, self.Int, self.errors
    

#    def plotT1rho()