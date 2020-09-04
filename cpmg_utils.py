# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 10:58:54 2020

Version 0.2 

- added automatic calculation of wdw length and placement
- moved the sweep of the window function into the cpmg_data class
- changed the window function inputs to meaningful values

@author: Selter
"""

import numpy as np
import nmrglue as ng
from scipy.fftpack import fft
from scipy.fftpack import fftshift
import scipy.signal
import matplotlib.pyplot as plt

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
        # print('cpmg object created')
        self.dpath = path+'/'+str(expno)
        self.dict, self.data = ng.bruker.read(self.dpath)
        self.SWH = self.dict['acqus']['SW_h']
        self.DW = 1.0/(2.0*self.SWH)
        self.TD = len(self.data)
        self.AQ = self.TD*2*self.DW
        self.SI = self.TD*4
        self.TDeff = self.TD
        self.time_axis = np.arange(0,self.AQ,2*self.DW)
        self.grpdly = int(round(self.dict['acqus']['GRPDLY']))+1
        self.n_echoes = self.dict['acqus']['L'][22]
        self.d3=self.dict['acqus']['D'][3]
        self.d6=self.dict['acqus']['D'][6]
        self.p1=self.dict['acqus']['P'][1]/1e6
        self.np_d6 = int(self.d6/(2*self.DW))
        self.np_echo = int((self.d3*2+self.d6+self.p1+2e-6)/(self.DW*2.0))
        self.windows = windows()
        print(self.np_echo)
        
    def set_SI(self,SI):
        self.SI = SI
        
    def set_TDeff(self,TDeff):
        self.TDeff = TDeff
        
    def define_sweep(self,plength,swh_pulse):
        self.sweep_length = plength
        self.sweep_width = swh_pulse
        
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

    def apply_em(self,lb):
        """
        Applies exponential linebroadening. The equation is slightly altered
        from the version used in Bruker Topspin to account for the difference
        in data point handling (real & imag alternating vs. complex data points)

        Parameters
        ----------
        lb : FLOAT or INTEGER
            LB constant analog to LB in Topspin.

        Returns
        -------
        None.

        """
        # print('EM linebroadeing applied')
        
        for n in range(len(self.data)):
            self.data[n] = self.data[n]*np.exp(-((n-1)*lb*np.pi)/(self.SWH))
        
        
    def apply_gm(self,gb,lb):
        """
        Applies gaussian linebroadening. This is similar to the 'gm' 

        Parameters
        ----------
        gb : FLOAT
            Value between 0 and 1, shift of the gaussian bell
            
        lb : FLOAT
            value describing the width of the bell curve, higher means narrower curve in the time domain.

        Returns
        -------
        None.

        """
        
        # self.T2 = 1/(lb)
        # self.T2np = (self.T2/(2*np.sqrt(2*np.log(2))))/(self.DW*2)
        # lb2 =  self.T2np 
        self.funky = getattr(windows,'gaussian')
        # self.funky2 = gauss(self.TD,lb2,self.TD*gb)
        self.data = self.funky(self.TD,lb,self.TD*gb)*self.data
                 

    def return_wdw_sweep(self,lb=10,start=0,end=1.0,res=1,**kwargs):
        
        wdw_type = kwargs.get('wdw','gaussian')
        self.SI = kwargs.get('SI', self.SI)
        self.TDeff = kwargs.get('TDeff', self.TDeff)
        apod = kwargs.get('apod',False )
        apod_type = kwargs.get('apod_type','em')
        apod_params = kwargs.get('apod_params',(100,0.1))
        
        self.sum = np.zeros(self.SI)
        self.start = int(self.np_d6*start)
        self.end = int(self.np_d6*end)
        self.range = int(self.np_d6*end)-int(self.np_d6*start)
        
        
        # This is the loop for shifting through the echo:
        for k in range(self.start,self.end,res):
            self.temp_gb = k/self.np_d6    
            
            
            self.apply_rpt_wdw(GB=self.temp_gb,LB=lb, n_echoes=self.n_echoes,wdw_name=wdw_type)    
                        
            if apod == True:
                if apod_type == 'em':
                    self.apply_em(lb=apod_params[0])
                if apod_type == 'gm':
                    self.apply_gm(lb=apod_params[0],gb=apod_params[1])
                    
            self.freq, self.fftdata = self.return_fft()
            self.sum = np.fmax(np.abs(self.sum),np.abs(self.fftdata))
            self.reload_fid()
            
            
        return self.freq,self.sum     
    
    
    def apply_rpt_wdw(self,**kwargs):
                      
        wdw_name = kwargs.get('wdw_name', 'gaussian')
        gb = kwargs.get('GB', 0.5)
        lb = kwargs.get('LB', 1000)
        n_echoes = kwargs.get('N_echoes',self.n_echoes)
        
        # LB is sigma in this case, which means lb (FWHM) is 2*LB
        
        self.center = self.grpdly + int(self.np_d6*gb)
        # self.T2 = 1/(lb)
        # self.T2np = (self.T2/(2*np.sqrt(2*np.log(2))))/(self.DW*2)
        # lb2 =  self.T2np
        self.funky = getattr(windows,wdw_name)
        self.wd = self.funky(npoints=self.np_echo,LB=lb,center=self.center)
        self.funky2 = getattr(windows,'rpt_wdw')
        self.wd2 = self.funky2(0,self.wd,n_echoes,0)
        self.apply_wdw(self.wd2)
            
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
    
    def return_wdw(self):
        return self.wd2
    
         
    def wdw_shift(self,shift):
        
        step = self.SWH/self.SI
        steps = int(shift*step)
        self.shifted_fftwdw = np.roll(self.fftwdw,steps)
        return self.shifted_fftwdw


###################
# Work in progress below
# 
##


    def new_proc(self,lb=10,start=0,end=1.0,res=1,**kwargs):
        
        wdw_type = kwargs.get('wdw','gaussian')
        self.SI = kwargs.get('SI', self.SI)
        self.TDeff = kwargs.get('TDeff', self.TDeff)
        apod = kwargs.get('apod',False )
        apod_type = kwargs.get('apod_type','em')
        apod_params = kwargs.get('apod_params',(100,0.1))
        
        self.sum = np.zeros(self.SI)
        self.start = int(self.np_d6*start)
        self.end = int(self.np_d6*end)
        self.range = int(self.np_d6*end)-int(self.np_d6*start)
        
        
        # This is the loop for shifting through the echo:
        for k in range(self.start,self.end,res):
            self.temp_gb = k/self.np_d6    
            
            self.current_window = self.apply
            self.fft_window = self.return_wdw_fft(self.wdwwdw, kwargs)
            
            self.apply_rpt_wdw(GB=self.temp_gb,LB=lb, n_echoes=self.n_echoes,wdw_name=wdw_type)    
                        
            if apod == True:
                if apod_type == 'em':
                    self.apply_em(lb=apod_params[0])
                if apod_type == 'gm':
                    self.apply_gm(lb=apod_params[0],gb=apod_params[1])
                    
            self.freq, self.fftdata = self.return_fft()
            self.sum = np.fmax(np.abs(self.sum),np.abs(self.fftdata))
            self.reload_fid()
            
            
        return self.freq,self.sum     

    def apply_rpt_wdw(self,**kwargs):
                      
        wdw_name = kwargs.get('wdw_name', 'gaussian')
        gb = kwargs.get('GB', 0.5)
        lb = kwargs.get('LB', 1000)
        n_echoes = kwargs.get('N_echoes',self.n_echoes)
        
        # LB is sigma in this case, which means lb (FWHM) is 2*LB
        
        self.center = self.grpdly + int(self.np_d6*gb)
        # self.T2 = 1/(lb)
        # self.T2np = (self.T2/(2*np.sqrt(2*np.log(2))))/(self.DW*2)
        # lb2 =  self.T2np
        self.funky = getattr(windows,wdw_name)
        self.wd = self.funky(npoints=self.np_echo,LB=lb,center=self.center)
        self.funky2 = getattr(windows,'rpt_wdw')
        self.wd2 = self.funky2(0,self.wd,n_echoes,0)
        self.apply_wdw(self.wd2)
      

      
    # def apply_wdw(self,wdw):
    #     """
    #     Multiplies each point of the FID with the value given in WDW

    #     Parameters
    #     ----------
    #     wdw : List of floats
    #         window function, must be as long as data set.

    #     Returns
    #     -------
    #     None.

    #     """
        
    #     self.rest = np.zeros(len(self.data)-len(wdw))
    #     self.wdw = np.concatenate((wdw, self.rest))
    #     self.data=self.data*self.wdw

# 
#  
################
              
         
    def return_fid(self):
        """
        return the time axis and the complex data points

        Returns
        -------
        time_axis: array of floats
            array of floats, corresponding to the time axis
        data: array of complex points
            data points, complex values

        """
               
        return self.time_axis, self.data
    
    def return_fft(self):
        """
        Returns the FFT of the raw data in the current processing state
        
        Parameters
        ----------
        SI : INTEGER
            zero filling to be applied.
        TDeff : INTEGER, optional
            TDeff to be applied. The default is 0.

        Returns
        -------
        ARRAY of FLOAT
            frequency axis.
        ARRAY of COMPLEX
            complex data points in the frequency domain.

        """
        
        if self.TDeff == 0:
            self.fftdata = fft(self.data, self.SI)
            self.fftdata = fftshift(self.fftdata)
            self.freq_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,self.SI)
            return self.freq_axis, self.fftdata        
        else:
            self.fftdata = fft(self.data[0:self.TDeff], self.SI)
            self.fftdata = fftshift(self.fftdata)
            self.freq_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,self.SI)
            return self.freq_axis, self.fftdata
    
    def return_wdw_fft(self,wdw,**kwargs):
        
        SI = kwargs.get('SI', self.SI)
        TDeff = kwargs.get('TDeff', 0)
        
        if self.TDeff == 0:
            self.fftwdw = fft(wdw, SI)
            self.fftwdw = fftshift(self.fftwdw)
            self.freq_wdw_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,SI)
            return self.freq_axis, self.fftwdw        
        else:
            self.fftwdw = fft(wdw[0:TDeff], SI)
            self.fftwdw= fftshift(self.fftwdw)
            self.freq_axis = np.linspace(-0.5*self.SWH,0.5*self.SWH,SI)
            return self.freq_axis, self.fftwdw
    
    
    
###############################################################################
# Definition of various window functions
# These are short factions called by the 
# CPMG_DATA class, but can also be used independently

# similarly, these are designed to be expandable
# and called by function name rather than hardcoded into the CPMG_DATA class.
# Therefore no changes need to be made to the class to make a new window function 
# available in the CPMG_DATA class.
#####

class windows(object):
    
    
    def __init__(self):
        """
        

        Returns
        -------
        None.

        """
    
    
    def gaussian(npoints=100,LB=10,center=50):
        """  
        Creates a Gaussian window function
    
        Parameters
        ----------
        npoints : INTEGER, optional
            number of points in the window. The default is 100.
        sigma : INTEGER, optional
            FWHH in points. The default is 10.
        center : INTEGER, optional
            where the center of the bell is. The default is 50.
    
        Returns
        -------
        out : list of floats
            returns a window function, in this case the gaussian bell.
    
        """
        lb2 = LB/(2*np.sqrt(2*np.log(2)))
        points = np.arange(npoints)
        out = np.exp(-0.5*((points-center)/(lb2))**2)
        return out          
            
    def lorentzian(npoints=100,width=10,center=50):
        """
        Creates the lorentzian window function
    
        Parameters
        ----------
        npoints : INTEGER, optional
            DESCRIPTION. The default is 100.
        width : INTEGER, optional
            DESCRIPTION. The default is 10.
        center : INTEGER, optional
            DESCRIPTION. The default is 50.
    
        Returns
        -------
        out : TYPE
            DESCRIPTION.
    
        """
        points = np.arange(npoints)
        x = (center-points)/(width/2)
        out = 1/(1+x**2)
        return out
    
    def rectangle(npoints=100,width=10,center=50):
        """
        
    
        Parameters
        ----------
        npoints : Integer, optional
            Length of total window function. The default is 100.
        width : INTEGER, optional
            width of the rectangle in points, if smaller than npoints zero padding is applied . The default is 10.
        center : INTEGER, optional
            Center position of the rectangle in points. The default is 50.
    
        Returns
        -------
        out : ARRAY of FLOAT
            DESCRIPTION.
    
        """
        out = np.zeros(npoints)
        for n in range(int(center-width*0.5),int(center+width*0.5)+1,1):
            out[n]=1.0
        
        return out
    
    
    def exponential(npoints=100,decay=0.01,center=50):
        """
        
    
        Parameters
        ----------
        npoints : Integer, optional
            Length of total window function. The default is 100.
        decay : INTEGER, optional
            decay of the Exponential in points, if smaller than npoints zero padding is applied . The default is 10.
        center : INTEGER, optional
            Center position of the rectangle in points. The default is 50.
    
        Returns
        -------
        out : ARRAY of FLOAT
            DESCRIPTION.
    
        """
        
        out = np.zeros(npoints)
        dec = np.zeros(npoints)
        decay2 = 1/decay
        for n in range(0,npoints):  
            dec[n]=1.0*np.exp(-decay2*n)
        # print(dec)
        out[center:] = dec[:(npoints-(center))]    
        flipped = np.flip(dec)
        out[:center+1] = flipped[(npoints-(center+1)):]
        
        # print(out)
        return out
    
    def tukey(npoints=100,width=10,center=50,alpha=0.5):
        """
        
    
        Parameters
        ----------
        npoints : Integer, optional
            Length of total window function. The default is 100.
        width : INTEGER, optional
            width of the function in points, if smaller than npoints zero padding is applied . The default is 10.
        center : INTEGER, optional
            Center position of the function in points. The default is 50.
    
        alpha : FLOAT, optional
            Range from 0.0 to 1.0. The default is 0.5.
    
        Returns
        -------
        out : ARRAY of FLOAT
            DESCRIPTION.
    
        """
        if (width % 2 == 0): 
            width=width
            window = scipy.signal.tukey(width,alpha,sym=False)
        # print(len(window))
            front = np.zeros(center-int(width*0.5))
        # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
        # print(len(end))
            out= np.concatenate((front, window, end))
        else: 
            width=width
            window = scipy.signal.tukey(width,alpha,sym=True)
            
        # print(len(window))
            front = np.zeros(center-int(width*0.5))
        # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
        # print(len(end))
            out= np.concatenate((front, window, end))
            
        return out
    
    def cheby(npoints=100,width=30,center=50):
        """
        
    
        Parameters
        ----------
        npoints : Integer, optional
            Length of total window function. The default is 100.
        width : INTEGER, optional
            width of the function in points, if smaller than npoints zero padding is applied . The default is 10.
        center : INTEGER, optional
            Center position of the function in points. The default is 50.
    
        Returns
        -------
        -------
        out : ARRAY of FLOAT
            DESCRIPTION.
        """
        if (width % 2 == 0):
            width= width
            window = scipy.signal.chebwin(width,100,sym=False)
            # print(len(window))
            front = np.zeros(center-int(width*0.5))
            # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
            # print(len(end))
            out= np.concatenate((front, window, end))
        else:
            width= width
            window = scipy.signal.chebwin(width,100,sym=True)
            # print(len(window))
            front = np.zeros(center-int(width*0.5))
            # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
            # print(len(end))
            out= np.concatenate((front, window, end))
        return out
    
    
    def flattop(npoints=100,width=30,center=50):
        """
        
    
        Parameters
        ----------
        npoints : Integer, optional
            Length of total window function. The default is 100.
        width : INTEGER, optional
            width of the function in points, if smaller than npoints zero padding is applied . The default is 10.
        center : INTEGER, optional
            Center position of the function in points. The default is 50.
    
        Returns
        -------
        -------
        out : ARRAY of FLOAT
            DESCRIPTION.
        """    
        if (width % 2 ==0):
            window = scipy.signal.flattop(width,sym=False)
        # print(len(window))
            front = np.zeros(center-int(width*0.5))
        # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
        # print(len(end))
            out= np.concatenate((front, window, end))
        else:
            window = scipy.signal.flattop(width,sym=False)
        # print(len(window))
            front = np.zeros(center-int(width*0.5))
        # print(len(front))
            end = np.zeros(npoints-int(width)-len(front))
        # print(len(end))
            out= np.concatenate((front, window, end))
        return out
    
    def sinc(npoints=100,lb=10,center=50):
        """
        Creates the sinc window function
        
        Parameters
        ----------
        npoints : INTEGER, optional
            number of points in the window. The default is 100.
        freq : INTEGER, optional
            frequency of the sine, higher values correspond to narrower window. The default is 10.
        center : INTEGER, optional
            center point of the window. The default is 50.
    
        Returns
        -------
        out : list of floats
            returns a window function,
    
        """
        freq=lb/(np.pi*10)
        print(freq)
        points = np.arange(npoints)
        out = np.zeros(npoints)
        for n in points:
            if n == center:
                out[n] = 1
            else:
                out[n]= np.sin(freq*(n-center))/(freq*(n-center))
        return out
        
    
    
    
    def rpt_wdw(front,wdw,nwdw,back):
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
