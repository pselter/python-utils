# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:02:58 2020

@author: philipp
"""


import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import scipy.signal
import math



##############################################################################
#
# Data fitting objects below
#
#
#



class TD_NMR_fit(object):
    
    def __init__(self,xData,yData):
        
        self.x = xData
        self.y = yData
        self.ncomp = 0
        self.nparams = 0
        # self.function = np.zeros(len(self.y))
        self.params = []
        self.func = np.zeros(len(self.y))
        self.fitdict = {
            "name": [], 
            "type": [],
            "parameter idx": [] ,
            "parameter name": [] ,
            "function": []
            }
        
    def gauss_comp(self,x,I0,T2):
        y = I0*np.exp(-0.5*(x/T2)**2)
        return y   
    
    def add_gauss(self,name,I0=1,T2=1e-3):

        self.fitdict["name"]+= [name]
        self.fitdict["type"] += ["gauss"]
        self.fitdict["parameter idx"] += [[self.nparams,self.nparams+1]]
        self.nparams += 2
        self.params += [I0,T2]
        self.fitdict["parameter name"] += [["Intensity","T2* constant"]]
        self.fitdict["function"] += [self.gauss_comp]
        self.ncomp += 1    

    def report_model(self):
        
        for n in range(self.ncomp):
            # print(self.fitdict)
            print("-----------------------------------------------")
            print("Component name : "+str(self.fitdict["name"][n]))
            print("component type : "+str(self.fitdict["type"][n]))
            print("parameter index : "+str(self.fitdict["parameter idx"][n]))
            l = 0
            for k in self.fitdict["parameter idx"][n]:
                print("Value of "+str(self.fitdict["parameter name"][n][l])+" : "+ str(self.params[k]))
                l += 1
            print("function" + str(self.fitdict["function"][n]))
            print("")
            # print("-----------------------------------------------")
            
    def function(self):
        
        for n in range(self.ncomp):
            paralist = []
            for k in self.fitdict["parameter idx"][n]:
                  paralist.append(self.params[k])
            print(paralist)
            self.func += self.fitdict["function"][n](self.x,*paralist)
        
        return self.func


##############################################################################
#
# Below is the data handling class
# includes loading raw data
# dealing with digital filter, group delay etc.
# 
# time shift for echos possible


class TD_NMR_data(object):
    
    def __init__(self,path,expno=1):
        """
        Creates the Time domain NMR data handling object

        Parameters
        ----------
        path : STRING
            path to the bruker data set.
        expno : INTEGER, optional
            EXPNO of the experiment. The default is 1.

        Returns
        -------
        None.

        """
        # print('cpmg object created')
        self.dpath = path+'/'+str(expno)
        self.dict, self.data = ng.bruker.read(self.dpath)
        self.SWH = self.dict['acqus']['SW_h']
        self.DW = 1.0/(2.0*self.SWH)
        self.TD = len(self.data)
        self.AQ = self.TD*2*self.DW
        self.time_axis = np.arange(0,self.AQ,2*self.DW)
        
    def return_fid(self):
        return self.time_axis, self.data
        
        
    def reload_fid(self):
        """
        Reloads the raw FID for processing
        Useful if something went wrong

        Returns
        -------
        None.

        """
        
        # print('loaded original FID')
        self.time_axis = np.arange(0,self.TD*2*self.DW,2*self.DW)
        self.dict, self.data = ng.bruker.read(self.dpath)
        
    def phase_fid(self,phase=90):
        self.phase = np.deg2rad(phase)
        self.data = self.data*np.exp(1j*self.phase)
        
    def convdta(self, verbose = True, shifted = False, add=0):
        self.reload_fid()
        
        # add plus one due to how python data structures work
        self.grpdly = int(round(self.dict['acqus']['GRPDLY']))+add+1
        
        if verbose == True:
            print("Group Delay points: "+str(self.grpdly))
        
        self.shifted_time_axis = self.time_axis[self.grpdly:] 
        
        if shifted == False:
            self.time_axis = self.time_axis[:-int(round(self.grpdly))]
        else:
            self.time_axis = self.shifted_time_axis
            
        self.data = self.data[int(round(self.grpdly)):]
        
        if verbose == True:
            print(len(self.time_axis))
            print(len(self.data))
    def truncate(self,TDeff):
        self.data = self.data[:TDeff]
        self.time_axis = self.time_axis[:TDeff]

        
    def correct_filter(self,shifted=False):
        self.reload_fid()
        
        self.grpdly = int(round(self.dict['acqus']['GRPDLY']))
        self.filt = self.data[:self.grpdly]
        self.zeros = np.zeros(len(self.filt))
        # print(self.filt)
        self.unfilt = self.zeros-self.filt
        # print(self.unfilt)
        self.unfilt = np.flip(self.unfilt)
        # print(self.unfilt)
        self.shifted_time_axis = self.time_axis[self.grpdly+1:]
        self.data = self.data[int(round(self.grpdly)+1):]        
        
        if shifted == False:
            self.time_axis = self.time_axis[:-int(round(self.grpdly+1))]
        else:
            self.time_axis = self.shifted_time_axis

        
        for n in range(len(self.unfilt)):
            self.data[n] = self.data[n]-self.unfilt[n]
            
        # return self.shifted_time_axis, self.data
        
    def timeshift(self,shift=0):
        self.time_axis = self.time_axis - shift
        
    def getnoise(self,mc=True,noisereg=1000):
        self.noise_avg = 0 
        if mc == True:
            for n in range(noisereg):
                self.noise_avg += self.mc_data[-n]    
                # print(self.mc_data[-n])
            self.noise_avg = self.noise_avg/noisereg
        else:
            self.noise_avg = np.std(self.data[-noisereg:])
        return self.noise_avg
  
    def calc_mag(self):
        self.mc_data = np.zeros(len(self.data))
        for n in range (len(self.data)):
            self.mc_data[n] = np.sqrt(self.data.real[n]**2+self.data.imag[n]**2)
        return self.time_axis, self.mc_data