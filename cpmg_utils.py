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

def create_sinc(npoints=100,width=10,center=50):
    
    points = np.arange(npoints)
    out = np.zeros(npoints)
    for n in points:
        if n == center:
            out[n] = 1
        else:
            out[n]= np.sin(width/100*(n-center))/(width/100*(n-center))
    return out
    

def create_gauss(npoints=100,width=10,center=50):
    """  
    Creates a Gaussian window function

    Parameters
    ----------
    npoints : INTEGER, optional
        number of points in the window. The default is 100.
    width : INTEGER, optional
        FWHH in points. The default is 10.
    center : INTEGER, optional
        where the center of the bell is. The default is 50.

    Returns
    -------
    out : list of floats
        returns a window function, in this case the gaussian bell.

    """
    points = np.arange(npoints)
    out = np.exp(-0.5*((points-center)/(width))**2)
    return out

def create_rp_wdw(front,wdw,nwdw,back):
    """
    Repeats a window function by concatenating it
    Optionally adds zeros front and back

    Parameters
    ----------
    front : INTEGER
        zeros to be added to the front.
    wdw : LIST of FLOATS
        previously created window function.
    nwdw : INTEGER
        number of wdw to be concatenated.
    back : INTEGER
        number of zeros to be added at the end.

    Returns
    -------
    out : LIST of FLOATS
        Concatenated list.

    """
    
    
    first = np.zeros(front)
    last = np.zeros(back)
    out = np.concatenate((first, wdw))
    for n in range(nwdw-1):
        out = np.concatenate((out, wdw))
    out = np.concatenate((out, last))
    return out

class cpmg_data(object):
    
    def __init__(self,path,expno=1):
        """
        Creates the cpmg_data handling object

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
        print('cpmg object created')
        self.dpath = path+'/'+str(expno)
        self.dict, self.data = ng.bruker.read(self.dpath)
        self.SWH = self.dict['acqus']['SW_h']
        self.DW = 1.0/(2.0*self.SWH)
        self.TD = len(self.data)
        self.AQ = self.TD*2*self.DW
        self.time_axis = np.arange(0,self.AQ,2*self.DW)
        
    def reload_fid(self):
        """
        Reloads the raw FID for processing
        Useful if something went wrong

        Returns
        -------
        None.

        """
        
        print('loaded original FID')
        self.time_axis = np.arange(0,self.TD*2*self.DW,2*self.DW)
        self.dict, self.data = ng.bruker.read(self.dpath)

    def apply_em(self,lb):
        """
        Applies exponential linebroadening

        Parameters
        ----------
        lb : FLOAT or INTEGER
            LB constant analog to LB in Topspin.

        Returns
        -------
        None.

        """
        print('EM linebroadeing applied')
        for n in range(len(self.data)):
            self.data[n] = self.data[n]*np.exp(-((n-1)*lb*np.pi)/(2.0*self.SWH))
        
    def apply_gm(self,gb,lb):
        """
         Applies gaussian linebroadening

        Parameters
        ----------
        gb : FLOAT
            Value between 0 and 1, shift of the gaussian bell
            
        lb : FLOAT
            negative value describing the width of the bell curve, higher means narrower curve.

        Returns
        -------
        None.

        """
        print('GM linebroadeing applied')
        self.b = lb/(2*gb*self.time_axis[-1])        
        self.data = self.data*np.exp((-1.0*lb*self.time_axis)-(-1.0*self.b*self.time_axis**2))
    
    def apply_wdw(self,wdw):
        """
        Multiplies each point of the FID with the value given in WDW

        Parameters
        ----------
        wdw : List of floats
            window function, must be as long as data set.

        Returns
        -------
        None.

        """
        
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
    