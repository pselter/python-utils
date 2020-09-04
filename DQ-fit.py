# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from matplotlib import pyplot as plt
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
def dual_expo_test(x,T1r_1,I1,T1r_2):
    I2=I1*2
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
###############################################################################   

class DQ_build_up(object):
    """General Object to handle DQ build up relaxation experiments
    
    expects the data being pre-processed via rowext in Topsin 
    or euqivalent treatment prior to use, also requires two datasets
        
    """
    
    
    #----------------------------------
    def __init__(self, path, expno_DQ,expno_REF,procno_DQ=100,procno_REF=100):
        """"Initialize the class
        
        nothing fancy, yet
        """
        
        self.path = path
        self.expno_DQ = expno_DQ
        self.expno_REF = expno_REF
        self.procno_DQ = procno_DQ
        self.procno_REF = procno_REF
        
    #----------------------------------
    
    
    def setup_exp(self,MAS=25000,increment=1,TD=32):
        """ 
        Procedure to set up the experimental parameters
        
        
        """
        self.MAS = MAS
        self.inc = increment
        self.TD = TD
        self.tau_rotor = 1/(self.MAS)
        self.time_axis = np.arange(self.tau_rotor,self.tau_rotor*self.TD,self.inc)
        
        return 
    
    
    
    #----------------------------------
    def getnoise(self,dat,noise_reg=(400,200)):
        """"Caculates the noise number based on the algorithm stated in the Bruker Topspin manual
            This is essentially the same as numpy.std(), which is why it will be replaced later on
            only retained for backwards compatibility
        """
        self.noiseidx0 = (np.abs(self.xaxis - noise_reg[0])).argmin()
        self.noiseidx1 = (np.abs(self.xaxis - noise_reg[1])).argmin()
        
        self.np_noise = abs(self.noiseidx1-self.noiseidx0)
        self.center = abs((self.noiseidx1-self.noiseidx0)*0.5)
        self.ynoise = sum(dat[self.noiseidx0:self.noiseidx1])
        self.y2noise = sum(dat[self.noiseidx0:self.noiseidx1]**2)
        self.xynoise = 0
              
        for n in range(1,int(self.center)):
            self.xynoise+=n*(dat[int(self.center+n)]-dat[int(self.center-n)])
        
        self.gnoise = np.sqrt((self.y2noise - (self.ynoise**2 + 3*self.xynoise**2/(self.np_noise**2-1))/self.np_noise)/(self.np_noise-1))
        print(self.gnoise)
        return self.gnoise
    
    
    def calcnoise(self,dat,noise_reg=(400,200)):
        """"Caculates the noise number based on  numpy.std(), which is essentially the same as the Bruker thing
        """
        self.nidx0 = (np.abs(self.xaxis - noise_reg[0])).argmin()
        self.nidx1 = (np.abs(self.xaxis - noise_reg[1])).argmin()
        self.noise = np.std(dat[self.nidx0:self.nidx1])
        return self.noise
    #----------------------------------



    def int_pseudo2d(self,region=(0,20),noise_reg=(400,200)):
        
        self.region = np.array(region)
        self.DQ_error = np.zeros((self.TD))
        self.REF_error = np.zeros((self.TD))
        self.integrals_DQ = np.zeros((self.TD))
        self.integrals_REF = np.zeros((self.TD))
        
        ###############################################################################
        
        #----------------------------------
        #Do the integration over all procnos
        for n in range(self.TD):
            
        #   load the data and find the spectral limits
        #----------------------------------
            self.spectrum_DQ = bruk.bruker1d(self.path,self.expno_DQ,procno=self.procno_DQ+n)
            self.spectrum_REF = bruk.bruker1d(self.path,self.expno_REF,procno=self.procno_REF+n)
            self.x_DQ, self.y_DQ = self.spectrum_DQ.plot1d()
            self.x_REF, self.y_REF = self.spectrum_REF.plot1d()
            self.xaxis_DQ = np.array(self.x_DQ)
            self.xaxis_REF = np.array(self.x_REF)
            self.zerofilling_factor_DQ,self.sn_fac_DQ = self.spectrum_DQ.calc_zfratio()
            self.zerofilling_factor_REF,self.sn_fac_REF = self.spectrum_REF.calc_zfratio()
                        
            #----------------------------------
            #check for the regions to be ordered correctly
            if self.region[0] < self.region[1]:
                tmp = self.region[0]
                self.region[0] = self.region[1]
                self.region[1] = tmp
                
            if noise_reg[0] < noise_reg[1]:
                tmp = noise_reg[0]
                noise_reg[0] = noise_reg[1]
                noise_reg[1] = tmp
            #---------------------------------- 
            
            #----------------------------------
            self.nidx0_DQ = (np.abs(self.xaxis_DQ - noise_reg[0])).argmin()
            self.nidx1_DQ = (np.abs(self.xaxis_DQ - noise_reg[1])).argmin()
            self.nidx0_REF = (np.abs(self.xaxis_REF - noise_reg[0])).argmin()
            self.nidx1_REF = (np.abs(self.xaxis_REF - noise_reg[1])).argmin()
            #----------------------------------
            self.noise_DQ = np.std(self.y_DQ[self.nidx0_DQ:self.nidx1_DQ])
            self.noise_REF = np.std(self.y_REF[self.nidx0_REF:self.nidx1_REF])
            
            #----------------------------------            
            #Do the integration for every defined region
            self.idx0_DQ = (np.abs(self.xaxis_DQ - self.region[0])).argmin()
            self.idx1_DQ = (np.abs(self.xaxis_DQ - self.region[1])).argmin()
            self.idx0_REF = (np.abs(self.xaxis_REF - self.region[0])).argmin()
            self.idx1_REF = (np.abs(self.xaxis_REF - self.region[1])).argmin()
            self.np_integral_DQ = self.idx1_DQ-self.idx0_DQ
            self.np_integral_REF = self.idx1_REF-self.idx0_REF
            #----------------------------------            

            #----------------------------------                        
            for p in range(self.idx1_DQ-self.idx0_DQ):
                self.integrals_DQ[n] = self.integrals_DQ[n] + self.y_DQ[self.idx0_DQ+p]/self.zerofilling_factor_DQ
            self.DQ_error[n] = self.noise_DQ*np.sqrt(self.np_integral_DQ)*self.sn_fac_DQ


            for p in range(self.idx1_REF-self.idx0_REF):
                self.integrals_REF[n] = self.integrals_REF[n] + self.y_REF[self.idx0_REF+p]/self.zerofilling_factor_REF
            self.REF_error[n] = self.noise_REF*np.sqrt(self.np_integral_REF)*self.sn_fac_REF
            #----------------------------------            

        




        #----------------------------------   
        ###############################################################################
        else:
            print(self.sl_times)
            return self.sl_times, self.integrals, self.sl_error




    def fit_pseudo2d(self,start_procno,peaklist,model,list_type='vd',fit_reg=(20,-10),use_sigma=False,noise_reg=(400,200),normalize=False,verbose=False,print_fits=False):
        
        self.sl_times = self.importvdlist(list_type)
        self.n_peaks = len(peaklist)
        self.npoints = len(self.sl_times)
        self.sl_int = np.zeros((self.npoints,self.n_peaks))
        self.sl_error = np.zeros((self.npoints,self.n_peaks))
        
        #define the fitting range
                       
        self.peaks = peaklist
        self.n_peaks = len(peaklist)
        self.integrals = np.zeros((self.npoints,self.n_peaks))
        
        ###############################################################################
        
        #----------------------------------
        #Do the integration over all procnos
        for n in range(self.npoints):
            
        #   load the data and find the spectral limits
            self.spectrum = bruk.bruker1d(self.path,self.expno,procno=start_procno+n)
            self.x, self.y = self.spectrum.plot1d()
            self.xaxis = np.array(self.x)
            self.idx0 = (np.abs(self.xaxis - fit_reg[0])).argmin()
            self.idx1 = (np.abs(self.xaxis - fit_reg[1])).argmin()
           
            self.zerofilling_factor,self.sn_fac = self.spectrum.calc_zfratio()
            
        #    find the noise region and calculate the sigma of the noise
        #    sigma of noise is used in the eventual fit as the error of each point
            if use_sigma ==True:    
                self.nidx0 = (np.abs(self.xaxis - noise_reg[0])).argmin()
                self.nidx1 = (np.abs(self.xaxis - noise_reg[1])).argmin()
                self.noise = np.std(self.y[self.nidx0:self.nidx1])
                self.sigmas = np.zeros(self.idx1-self.idx0)
                for k in range(len(self.sigmas)):
                    self.sigmas[k]=self.noise
                
        #   Do the fit and save the results     
            self.initguess = np.zeros(self.n_peaks)
            self.initguess = self.initguess+100
            if use_sigma == True:
                self.peak_popt, self.peak_pcov = curve_fit(model,self.x[self.idx0:self.idx1],self.y[self.idx0:self.idx1],sigma=self.sigmas,absolute_sigma=True,p0=self.initguess,maxfev=5000)
            else:
                self.peak_popt, self.peak_pcov = curve_fit(model,self.x[self.idx0:self.idx1],self.y[self.idx0:self.idx1],p0=self.initguess,maxfev=5000)
            self.peak_fitted = model(self.x[self.idx0:self.idx1],*self.peak_popt)
            self.peak_residual = self.y[self.idx0:self.idx1]-self.peak_fitted
            self.sl_error[n,:] = np.sqrt(np.diagonal(self.peak_pcov))/self.peak_popt
            
        #   Do the integration of the individual peaks
        #   NOTE: the factors for zerofilling and sn_fac are empirical factors to account for differences in zerofilling
            for k in range(self.n_peaks):
                for p in range(self.idx1-self.idx0):
                    self.integrals[n,k] = self.integrals[n,k] + self.peaks[k](self.x[self.idx0+p],self.peak_popt[k])/self.zerofilling_factor
            self.sl_error[n,:] = self.sl_error[n,:]*self.integrals[n,:]*self.sn_fac
        
        
        #   Some optional output to better track what is happening    
            if verbose == True:
                print('---------------')
                print('--Integrals--')
                for k in range(self.n_peaks):
                    print(self.integrals[n,k])
                print('---------------')
        #   optional output of all the fits and spectra to track fitting
            if print_fits == True:
                fig = plt.figure()
                plt.xlim(fit_reg)
                if use_sigma == True:
                    plt.errorbar(self.x[self.idx0:self.idx1],self.y[self.idx0:self.idx1],yerr=self.sigmas,label='exp')
                else:
                    plt.errorbar(self.x[self.idx0:self.idx1],self.y[self.idx0:self.idx1],label='exp')
                plt.plot(self.x[self.idx0:self.idx1],self.peak_fitted,label='fit')
                plt.plot(self.x[self.idx0:self.idx1],self.peak_residual,label='diff')
                for k in range(self.n_peaks):
                    plt.plot(self.x[self.idx0:self.idx1],self.peaks[k](self.x[self.idx0:self.idx1],self.peak_popt[k]))
                plt.show()
                
        if normalize == True:
            self.norm_integrals = np.zeros_like(self.integrals)
            self.norm_sl_error = np.zeros_like(self.sl_error)
            for k in range(self.n_peaks):
                self.norm_integrals[:,k] = self.integrals[:,k]/self.integrals[:,k].max()*100
                self.norm_sl_error[:,k] = self.sl_error[:,k]/self.integrals[:,k].max()*100
            return self.sl_times, self.norm_integrals, self.norm_sl_error
        #----------------------------------   
        ###############################################################################
        else:
            print(self.sl_times)
            return self.sl_times, self.integrals, self.sl_error


###############################################################################
#    not sure if these actually work at this point

#    def T1rho2(self,start_procno,npoints,region=(20,-10),noise_reg=(400,200), normalize=True):
#        """Function to do the integration, returns the integral and time in an array         
#        written with T1_rho measurements by spin-locking in mind
#        
#        assumes a vplist being present
#        - start_procno :the first procno with slices
#        - npoints      :the number of procnos to evaluate
#        - region       :the region to integrate in ppm
#        - normalize    :normalizes the data to 100
#        """
#        
#        self.vplist = self.importvdlist(list_type='vp')
#        self.intensity = np.zeros(npoints)
#        self.noise = np.zeros(npoints)
#        self.sn = np.zeros(npoints)
#        self.data = np.zeros((npoints,3))
#
#        
#        #----------------------------------
#        #Do the integration over all procnos
#        for n in range(npoints):
#            spectrum = bruk.bruker1d(self.path,self.expno,procno=start_procno+n)
#            self.x, self.y = spectrum.plot1d()
#            self.xaxis = np.array(self.x)
#            
#            
#            self.idx0 = (np.abs(self.xaxis - region[0])).argmin()
#            self.idx1 = (np.abs(self.xaxis - region[1])).argmin()
#            self.intensity[n] = np.amax(self.y[self.idx0:self.idx1])
#            self.noise[n]=self.getnoise(self.y)
#            self.sn[n]=self.intensity[n]/(2*self.noise[n])
#        #----------------------------------
#        
#        #----------------------------------
#        # Check for normalization
#        if normalize == True:
#            self.norm_int = self.intensity/self.intensity.max()*100
#            self.data[:,0] = self.vplist
#            self.data[:,1] = self.norm_int
#            self.data[:,2] = 1/self.sn
#        else:
#            self.data[:,0] = self.vplist
#            self.data[:,1] = self.intensity
#            self.data[:,2] = 1/self.sn
#        #----------------------------------
#        return self.data
#
#
#    
#    def T1rho(self,start_procno,npoints,region=(20,-10), normalize=True):
#        """Function to do the integration, returns the integral and time in an array         
#        written with T1_rho measurements by spin-locking in mind
#        
#        assumes a vplist being present
#        - start_procno :the first procno with slices
#        - npoints      :the number of procnos to evaluate
#        - region       :the region to integrate in ppm
#        - normalize    :normalizes the data to 100
#        """
#        
#        self.vplist =  self.importvdlist(list_type='vp')
#        self.integral = np.zeros(npoints)
#        #self.time = np.arange(0,self.vplist[-1],self.vplist[0])
#        self.data = np.zeros((npoints,2))
#        #self.res = np.zeros((npoints,2))
#        #self.fit =  np.zeros((len(self.time),2))
#        
#        #----------------------------------
#        #Do the integration over all procnos
#        for n in range(npoints):
#            spectrum = bruk.bruker1d(self.path,self.expno,procno=start_procno+n)
#            self.x, self.y = spectrum.plot1d()
#            self.xaxis = np.array(self.x)
#            self.idx0 = (np.abs(self.xaxis - region[0])).argmin()
#            self.idx1 = (np.abs(self.xaxis - region[1])).argmin()
#            for l in range(self.idx1-self.idx0):
#                self.integral[n] += self.y[self.idx0+l]
#        #----------------------------------
#        
#        #----------------------------------
#        # Check for normalization
#        if normalize == True:
#            self.norm_int = self.integral/self.integral.max()*100
#            print(self.vplist)
#            self.data[:,0] = self.vplist
#            self.data[:,1] = self.norm_int
#        else:
#            self.data[:,0] = self.vplist
#            self.data[:,1] = self.integral
#        #----------------------------------
#        return self.data


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
            
        if n_expo == 99:
            self.popt, self.pcov = curve_fit(dual_expo_test,self.data[:,0],self.data[:,1],maxfev=5000,p0=(1.0e-2,33,1.0e-3))
            self.fit[:,1] = dual_expo_test(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-dual_expo_test(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
#            self.Int[1] = self.popt[3]

        return self.fit, self.res, self.T2, self.Int, self.errors



    def exp_decay2(self,data,n_expo=1):
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
            self.popt, self.pcov = curve_fit(single_expo,self.data[:,0],self.data[:,1],sigma=self.data[:,2],maxfev=5000,p0=(1.0e-2,100))
            self.fit[:,1] = single_expo(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-single_expo(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2 = self.popt[0]
            self.Int = self.popt[1]
    
        if n_expo == 2:
            self.popt, self.pcov = curve_fit(dual_expo,self.data[:,0],self.data[:,1],sigma=self.data[:,2],maxfev=5000,p0=(1.0e-2,50,1.0e-3,50))
            self.fit[:,1] = dual_expo(self.time,*self.popt)
            self.res[:,1] = self.data[:,1]-dual_expo(self.data[:,0],*self.popt)
            
            self.errors = np.sqrt(np.diagonal(self.pcov))
            self.T2[0] = self.popt[0]
            self.Int[0] = self.popt[1]
            self.T2[1] = self.popt[2]
            self.Int[1] = self.popt[3]
        
        if n_expo == 3:
            self.popt, self.pcov = curve_fit(triple_expo,self.data[:,0],self.data[:,1],sigma=self.data[:,2],maxfev=5000,p0=(5.0e-2,33,1.0e-2,33,5.0e-3,33))
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