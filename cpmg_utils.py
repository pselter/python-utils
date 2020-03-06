# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:58:54 2020

@author: Selter
"""

import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import math



def create_gauss(npoints=100,width=10,center=50):
    points = np.arange(npoints)
    out = np.exp(-0.5*((points-center)/(width))**2)
    return out

def create_rp_wdw(front,wdw,nwdw,back):
    first = np.zeros(front)
    last = np.zeros(back)
    out = np.concatenate((first, wdw))
    for n in range(nwdw-1):
        out = np.concatenate((out, wdw))
    out = np.concatenate((out, last))
    return out

class cpmg_data(object):
    
    def __init__(self,path,expno=1):
        print('cpmg object created')
        self.dpath = path+'/'+str(expno)
        self.dict, self.data = ng.bruker.read(self.dpath)
        self.SWH = self.dict['acqus']['SW_h']
        self.DW = 1.0/(2.0*self.SWH)
        self.TD = len(self.data)
        self.AQ = self.TD*2*self.DW
        self.time_axis = np.arange(0,self.AQ,2*self.DW)
        
    def reload_fid(self):
        print('loaded original FID')
        self.time_axis = np.arange(0,self.TD*2*self.DW,2*self.DW)
        self.dict, self.data = ng.bruker.read(self.dpath)

    def apply_em(self,lb):
        print('EM linebroadeing applied')
        for n in range(len(self.data)):
            self.data[n] = self.data[n]*np.exp(-((n-1)*lb*np.pi)/(2.0*self.SWH))
        
    def apply_gm(self,gb,lb):
        print('GM linebroadeing applied')
        self.b = lb/(2*gb*self.time_axis[-1])        
        self.data = self.data*np.exp((-1.0*lb*self.time_axis)-(-1.0*self.b*self.time_axis**2))
    
    def apply_wdw(self,wdw):
        self.rest = np.zeros(len(self.data)-len(wdw))
        self.wdw = np.concatenate((wdw, self.rest))
        self.data=self.data*self.wdw
         
    def return_fid(self):
        return self.time_axis, self.data
    
    def return_fft(self,SI,TDeff=0):
        if TDeff == 0:
            self.fftdata = fft(self.data, SI)
            self.fftdata = fftshift(self.fftdata)
            self.freq_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,SI)
            return self.freq_axis, self.fftdata        
        else:
            self.fftdata = fft(self.data[0:TDeff], SI)
            self.fftdata = fftshift(self.fftdata)
            self.freq_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,SI)
            return self.freq_axis, self.fftdata
    