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
       
    def calc_mag(self):
        
        self.mc_data = np.zeros(len(self.data))
        for n in range (len(self.data)):
            self.mc_data[n] = np.sqrt(self.data.real[n]**2+self.data.imag[n]**2)
        return self.time_axis, self.mc_data