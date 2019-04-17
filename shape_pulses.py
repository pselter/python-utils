#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 14:06:49 2018

@author: akhansen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.fftpack import fftshift
import math


def plot_shape(time,amp,freq,phase_prog):
    """plot the pulse shape.
    
    time        -- time axis
    amp         -- amp list
    freq        -- freq list
    phase_prog  -- phase list
    """  
    fig = plt.figure(figsize=(10,7))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1, rowspan=1)
    plt.plot(time,freq)
    plt.title('Frequency')
    ax1.set_ylabel(r'Frequency / kHz', fontname="Arial", fontsize=12)
    ax1.set_xlabel(r'Time / µs', fontname="Arial", fontsize=12)
    
    ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
    plt.plot(time,amp)
    plt.title('Amplitude')
    ax2.set_ylabel(r'Amplitude / a.u.', fontname="Arial", fontsize=12)
    ax2.set_xlabel(r'Time / µs', fontname="Arial", fontsize=12)
    
    ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
    plt.plot(time,phase_prog)
    #plt.scatter(time,phase_prog)
    plt.title('Phase')
    ax3.set_ylabel(r'Phase / degree', fontname="Arial", fontsize=12)
    ax3.set_xlabel(r'Time / µs', fontname="Arial", fontsize=12)
    
    ax4 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)
    fid = np.zeros(int(len(time)),dtype=complex)
    for n in range(0,len(time)):
        fid[n] = amp[n]*np.exp(1j*phase_prog[n]*np.pi/180)
    for x in range(len(fid)):
        plt.plot([0,fid[x].real],[0,fid[x].imag],'ro',label='python')
    limit=np.max(np.ceil(np.absolute(fid))) # set limits for axis
    plt.xlim((-limit,limit))
    plt.ylim((-limit,limit))
    plt.title('Phase')
    ax4.set_ylabel(r'Imaginary', fontname="Arial", fontsize=12)
    ax4.set_xlabel(r'Real', fontname="Arial", fontsize=12)
    
    plt.show()
    plt.close()
    return;

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx];

    
def analysize_shape(time,amp,freq,phase_prog): 
    """create the actual pulse and FFT.
    
    time        -- time axis
    amp         -- amp list
    freq        -- freq list
    phase_prog  -- phase list
    """
    profile = np.zeros((int(len(time)),2))
    fid = np.zeros(int(len(time)),dtype=complex)
    for n in range(0,len(profile)):
        profile[n,0] = amp[n]*np.cos(phase_prog[n]*np.pi/180)
        profile[n,1] = amp[n]*np.sin(phase_prog[n]*np.pi/180)
    for n in range(0,len(fid)):
        fid[n] = amp[n]*np.exp(1j*phase_prog[n]*np.pi/180)
        
    fig = plt.figure(figsize=(10,4))
    plt.subplots_adjust(wspace=0.4, hspace=0.2)
    ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=1, rowspan=1)
    plt.plot(time,fid.real,label='Real')
    plt.plot(time,fid.imag,label='Imag')
    plt.title('Pulse')
    plt.legend()
    ax1.set_ylabel(r'Intensity / a.u.', fontname="Arial", fontsize=12)
    ax1.set_xlabel(r'Time / µs', fontname="Arial", fontsize=12)
    
    ax2 = plt.subplot2grid((1, 2), (0, 1), colspan=1, rowspan=1)
    plt.title('Fourier Transform')
    SI=4*int(len(time))
    step_size = time[-1]/int(len(time))
    fx= np.linspace(-0.5*(1/(step_size)),0.5*(1/(step_size)),SI)
    fy = fft(fid,n=SI)
    fys = fftshift(fy)
    
    fid_abs = np.zeros(int(SI))
    for n in range(0,len(fid_abs)):
        fid_abs[n] = np.sqrt(fys[n].real**2+fys[n].imag**2)
        
        
    maximum = np.amax(fid_abs) 
    
    max_half = np.zeros(int(SI))
    for n in range(0,len(max_half)):
        max_half[n] = maximum/2
        
    idx = np.argwhere(np.diff(np.sign(max_half - fid_abs)) != 0).reshape(-1) + 0
    FWHH=fx[idx[1]]-fx[idx[0]]
    print('FWHH: '+str(np.round(FWHH*1000))+' kHz')
    
    plt.plot(fx,fys.real, label='Real')
    plt.plot(fx,fys.imag, label='Imag')
    plt.plot(fx,fid_abs, label='Abs')
#    plt.plot((-1,1),(maximum/2,maximum/2), label='HH')
    plt.plot(fx[idx], max_half[idx], 'ro')
    plt.xlim((-1,1))
    plt.legend()
    ax2.set_ylabel(r'Intensity /a.u.', fontname="Arial", fontsize=12)
    ax2.set_xlabel(r'Frequency', fontname="Arial", fontsize=12)
    
    plt.show()
    plt.close()
    return;

def calc_SHAP_par(delta_sigma,MAS,tan_k=20,c=10):
    """calculate the requirements for a SHAP pulse.
    
    delta_sigma -- maximum offset in kHz
    MAS         -- MAS freq.in kHz
    tan_k       -- shape parameter
    c           -- shape parameter
    """
    
    d_sigma =delta_sigma*2*np.pi
    wMAS = MAS*2*np.pi
    w_max_g = np.sqrt(1.41*d_sigma*wMAS)
    k = np.arctan(tan_k)
    print('Nutation frequency must be larger than: '+str(w_max_g/(2*np.pi)))
    delta_w_soll = 0.5*tan_k*d_sigma
    print('maximum offset at least: '+str( delta_w_soll/(2*np.pi)))
    tp_soll = ((2*delta_w_soll)/(d_sigma*wMAS))*(k/tan_k)
    print('pulse length at least '+str( tp_soll*1000))
    return;

def create_SHAP(pulse_length,step_size,sweep_width,tan_k=20,c=10,create=False):
    """create a SHAP pulse shape and returns time,amp,freq,phase lists.
    
    pulse_length -- pulse length in µs
    step_size    -- step_size in µs (e.g. 0.1 or 0.05)
    sweep_width  -- sweep_width in kHz
    tan_k        -- shape parameter (default 20)
    c            -- shape parameter (default 10)
    create       -- write shape file? (default False)
    """
    
    ##length in µs 
    tp=pulse_length

    ## stepsize in µs
    step_size= step_size

    ## maximum offset in kHz
    delta_w = sweep_width/2

    ## shape parameters:
    tan_k = tan_k
    k = np.arctan(tan_k)
    print(k)
    c = c


    ##calculate number of points in shape
    nump= int(tp/step_size)
    ## maximum amplitude
    w_max = 100
    ##create time axes
    t_1 = np.linspace(0,tp*0.5,int(nump/2))
    t_2 = np.linspace(tp*0.5,tp,int(nump/2)+1)
    t = np.linspace(0.0,tp,nump+1)
    ### create amplitude profile
    w1_1 =  w_max*np.tanh(2*c*t_1/tp)
    w1_2 =  w_max*np.tanh(2*c*(1-(t_2/tp)))
    w1 = np.concatenate((w1_1,w1_2))
    ## calculate frequency profile
    d_w = delta_w*(np.tan(k*(1-(2*t/tp)))/np.tan(k)) 
    ## calculate the corresponding phase profile
    phase = np.zeros(int(nump+1))
    ## uses analytical expression
    phase = (tp*delta_w*(1/tan_k)*np.log(np.cos(k-(2*k*t)/tp)))/(2*k)
    
    ### old buggy versions
#    for n in range(0,len(phase)):
#        phase[n] = d_w[n]*((t[n]-tp/2)/1000)*360
#        phase[n] = phase[n-1]+(d_w[n]*1000)*(step_size/1e6)*360
 
    
    ### wrap phase
    phase_prog = np.zeros(int(nump+1))
    for n in range(0,len(phase)):
            phase_prog[n] = math.fmod(phase[n],360)  
    
    ####writing shape file
    if create == True:
        filename = 'SHAP_'+str(tp)+'us_'+str(delta_w*2)+'_'+str(step_size)
    
        outfile = open(filename,'w')
        outfile.write('##TITLE= SHAP-'+str(tan_k)+' tan_k (duration '+str(tp)+' us, sweep width '+str(delta_w*2)+' kHz, step size: '+str(step_size)+' us'+'\n')
        outfile.write('##USAGE= SHAP pulse test'+'\n')
        outfile.write('##JCAMP-DX= 5.00 $$ Bruker JCAMP library \n##DATA TYPE= Shape Data \n##ORIGIN= Generated from wurst program \n##DATE= 02/22/18  \n##TIME= 16:02:37 \n')
        outfile.write('##MINX= 0.000000e+00 \n##MAXX= 1.000000e+02 \n##MINY= 0.000000e+00 \n##MAXY= 3.600000e+02 \n##$SHAPE_EXMODE= Adiabatic \n##$SHAPE_TOTROT= 1.800000e+02 \n##$SHAPE_TYPE= Inversion \n##$SHAPE_BWFAC= 0.000000e+00 \n##$SHAPE_INTEGFAC= 0.000000e+00 \n##$SHAPE_MODE= 1 \n')
        outfile.write('##NPOINTS= '+str(nump+1)+'\n'+'##XYPOINTS= (XY..XY) \n')
        for n in range(0,len(phase_prog)):
            outfile.write(str(w1[n])+' '+str(phase_prog[n])+'\n')
        outfile.write('##END=') 
        outfile.close()
        
        print('shape file for a SHAP pulse was created\n')
        print('Filename: '+str(filename)+'\n')
        print('Parameters:\n')
        print('Pulse length: '+str(pulse_length)+'\n')
        print('Sweep width: '+str(sweep_width)+'\n')
        print('number of points: '+str(nump+1)+'\n')
        print('step_size: '+str(step_size)+'\n')
    return t,w1,d_w,phase_prog;



def create_BRATWURST(pulse_length,step_size,sweep_width,N=80,create=False):
    """create a BRATWURST pulse shape file and returns time,amp,freq,phase lists.
    
    pulse_length -- pulse length in µs
    step_size    -- step_size in µs (e.g. 0.1 or 0.05)
    sweep_width  -- sweep_width in kHz
    N            -- shape parameter (default 80)
    create       -- write shape file? (default False)
    """

    ## units kHz
    sweep_width = sweep_width
    ## pulse in µs 
    pulse_length = pulse_length
    ##step size in µs
    t_inc = step_size
    nump = pulse_length/t_inc
    ## pulse power in kHz
    b1_max = 100
    N = N
    
    nominal_sweep_rate  = sweep_width/pulse_length
    #print('nominal sweep rate: '+str(nominal_sweep_rate)+' kHz/µs')
    time = np.linspace(0,pulse_length,int(nump+1))
    #print(time)
    
    offset_start = (-1)*sweep_width*0.5
    offset_end = offset_start+sweep_width
      
    freq = offset_start+(nominal_sweep_rate*time)
    amp = b1_max*(1-np.abs(np.cos((np.pi*time)/pulse_length))**N)
    q = np.abs(amp**2/(nominal_sweep_rate*1e9))
    
    current_sweep_rate = np.zeros(int(nump+1))
    for n in range(0,len(current_sweep_rate)-1):
        current_sweep_rate[n] = (freq[n+1]-freq[n])/t_inc
    current_sweep_rate[-1] = current_sweep_rate[0]
       
    sweep_corr_factor = (amp**2)/(b1_max**2)
    corr_sweep_rate = current_sweep_rate*sweep_corr_factor
    q_current = np.abs(amp**2/(2*np.pi*current_sweep_rate*1e9))
    corr_freq = np.zeros(int(nump+1))
    for n in range(0,len(corr_freq)):
        corr_freq[n] = corr_freq[n-1]+corr_sweep_rate[n]*t_inc
    zero_off = corr_freq[-1]*0.5
    for n in range(0,len(corr_freq)):
        corr_freq[n] = corr_freq[n]-zero_off
    phase = np.zeros(int(nump+1))
    for n in range(0,len(phase)):
        phase[n] = phase[n-1]+(corr_freq[n]*1000)*(t_inc/1e6)*360
#        phase[n] = (corr_freq[n])*((time[n]-pulse_length/2)/1e3)*180
    zero_phase = phase[int(nump/2)]
       
    for n in range(0,len(corr_freq)):
        phase[n] = phase[n]-zero_phase
        phase_prog = np.zeros(int(nump+1))
    for n in range(0,len(phase)):
        phase_prog[n] = math.fmod(phase[n],360)  
   
    q_current = np.abs(amp**2/(corr_sweep_rate))    
    
    ####writing shape file
    if create == True:
        filename = 'Bratwurst'    
        
        outfile = open(filename,'w')
        outfile.write('##TITLE= WURST-'+str(N)+' shape (duration '+str(pulse_length)+' us, sweep width '+str(sweep_width)+' kHz, step size: '+str(t_inc)+' us, sweep rate optimized'+'\n')
        outfile.write('##USAGE= WURST sweeps of the satellite transitions with symmetric offsets at 1300.0 kHz'+'\n')
        outfile.write('##JCAMP-DX= 5.00 $$ Bruker JCAMP library \n##DATA TYPE= Shape Data \n##ORIGIN= Generated from wurst program \n##DATE= 02/22/18  \n##TIME= 16:02:37 \n')
        outfile.write('##MINX= 0.000000e+00 \n##MAXX= 1.000000e+02 \n##MINY= 0.000000e+00 \n##MAXY= 3.600000e+02 \n##$SHAPE_EXMODE= Adiabatic \n##$SHAPE_TOTROT= 1.800000e+02 \n##$SHAPE_TYPE= Inversion \n##$SHAPE_BWFAC= 0.000000e+00 \n##$SHAPE_INTEGFAC= 0.000000e+00 \n##$SHAPE_MODE= 1 \n')
        outfile.write('##NPOINTS= '+str(nump+1)+'\n'+'##XYPOINTS= (XY..XY) \n')
        for n in range(0,len(corr_freq)):
            outfile.write(str(amp[n])+' '+str(phase_prog[n])+'\n')
        outfile.write('##END=') 
        outfile.close()
        
        print('shape file for a BRATWURST pulse was created\n')
        print('Filename: '+str(filename)+'\n')
        print('Parameters:\n')
        print('Pulse length: '+str(pulse_length)+'\n')
        print('Sweep width: '+str(sweep_width)+'\n')
        print('number of points: '+str(nump)+'\n')
        print('step_size: '+str(step_size)+'\n')
        print('N: '+str(N)+'\n')
    return time,amp,corr_freq,phase_prog;


    
    
def create_WURST(pulse_length,step_size,sweep_width,N=80,create=False):
    """create a WURST pulse shape file and returns time,amp,freq,phase lists.
    
    pulse_length -- pulse length in µs
    step_size    -- step_size in µs (e.g. 0.1 or 0.05)
    sweep_width  -- sweep_width in kHz
    N            -- shape parameter (default 80)
    create       -- write shape file? (default False)
    """

    ## units kHz
    sweep_width = sweep_width
    ## pulse in µs 
    pulse_length = pulse_length
    ##step size in µs
    t_inc = step_size
    nump = pulse_length/t_inc
    ## pulse power in kHz
    b1_max = 100
    N = N
    
    nominal_sweep_rate  = sweep_width/pulse_length
    #print('nominal sweep rate: '+str(nominal_sweep_rate)+' kHz/µs')
    time = np.linspace(0,pulse_length,int(nump))
    #print(time)
    
    offset_start = (-1)*sweep_width*0.5
    offset_end = offset_start+sweep_width
      
    freq = offset_start+(nominal_sweep_rate*time)
    amp = b1_max*(1-np.abs(np.cos((np.pi*time)/pulse_length))**N)
    q = np.abs(amp**2/(nominal_sweep_rate*1e9))
        
    phase = np.zeros(int(nump))
    for n in range(0,len(phase)):
        phase[n] = (freq[n])*((time[n]-pulse_length/2)/1e3)*180
    zero_phase = phase[int(nump/2)]
       
    for n in range(0,len(phase)):
        phase[n] = phase[n]-zero_phase
        
    phase_prog = np.zeros(int(nump))
    for n in range(0,len(phase)):
        phase_prog[n] = math.fmod(phase[n],360)  
       
    ####writing shape file
    if create == True:
        filename = 'WURST'    
        
        outfile = open(filename,'w')
        outfile.write('##TITLE= WURST-'+str(N)+' shape (duration '+str(pulse_length)+' us, sweep width '+str(sweep_width)+' kHz, step size: '+str(t_inc)+' us, sweep rate optimized'+'\n')
        outfile.write('##USAGE= WURST sweeps of the satellite transitions with symmetric offsets at 1300.0 kHz'+'\n')
        outfile.write('##JCAMP-DX= 5.00 $$ Bruker JCAMP library \n##DATA TYPE= Shape Data \n##ORIGIN= Generated from wurst program \n##DATE= 02/22/18  \n##TIME= 16:02:37 \n')
        outfile.write('##MINX= 0.000000e+00 \n##MAXX= 1.000000e+02 \n##MINY= 0.000000e+00 \n##MAXY= 3.600000e+02 \n##$SHAPE_EXMODE= Adiabatic \n##$SHAPE_TOTROT= 1.800000e+02 \n##$SHAPE_TYPE= Inversion \n##$SHAPE_BWFAC= 0.000000e+00 \n##$SHAPE_INTEGFAC= 0.000000e+00 \n##$SHAPE_MODE= 1 \n')
        outfile.write('##NPOINTS= '+str(nump+1)+'\n'+'##XYPOINTS= (XY..XY) \n')
        for n in range(0,len(freq)):
            outfile.write(str(amp[n])+' '+str(phase_prog[n])+'\n')
        outfile.write('##END=') 
        outfile.close()
        
        print('shape file for a BRATWURST pulse was created\n')
        print('Filename: '+str(filename)+'\n')
        print('Parameters:\n')
        print('Pulse length: '+str(pulse_length)+'\n')
        print('Sweep width: '+str(sweep_width)+'\n')
        print('number of points: '+str(nump)+'\n')
        print('step_size: '+str(step_size)+'\n')
        print('N: '+str(N)+'\n')
    return time,amp,freq,phase_prog;
  
    
def create_RECT(pulse_length,step_size,offset=0,ph=0,padding=True):
    """create a RECT pulse shape file and returns time,amp,freq,phase lists.
    
    pulse_length -- pulse length in µs
    step_size    -- step_size in µs (e.g. 0.1 or 0.05)
    offset       -- sweep_width in kHz (default 0)
    ph           -- phase in degree (default 0)
    padding      -- zeropadding before and aft (default True)
    """
    if padding == True:
        ## pulse in µs 
        pulse_length = pulse_length*4
        ##step size in µs
        t_inc = step_size
        nump = pulse_length/t_inc
        ## pulse power in kHz
        b1_max = 100
            
        nominal_sweep_rate  = 0
        #print('nominal sweep rate: '+str(nominal_sweep_rate)+' kHz/µs')
        time = np.linspace(0,pulse_length,int(nump))
        #print(time)
        
        freq = offset+time*0.0
        
        amp= np.zeros(int(nump))
        for n in range(0,int(3*nump/8)):
            amp[n] = time[n]*0.0
        for n in range(int(3*nump/8),int(5*nump/8)):
            amp[n] = b1_max+time[n]*0.0
        for n in range(int(5*nump/8),int(nump)):
            amp[n] = time[n]*0.0
            
        phase = np.zeros(int(nump))
        for n in range(0,len(phase)):
            phase[n] = (freq[n])*((time[n]-pulse_length/2)/1e3)*360
     
        phase_prog = np.zeros(int(nump))
        for n in range(0,len(phase)):
            phase_prog[n] = math.fmod(phase[n],360)  
            
            
    if padding == False:
        ## pulse in µs 
        pulse_length = pulse_length
        ##step size in µs
        t_inc = step_size
        nump = pulse_length/t_inc
        ## pulse power in kHz
        b1_max = 100
            
        nominal_sweep_rate  = 0
        #print('nominal sweep rate: '+str(nominal_sweep_rate)+' kHz/µs')
        time = np.linspace(0,pulse_length,int(nump))
        #print(time)
        
        
          
        freq = offset+time*0.0
        phase = np.zeros(int(nump))
        for n in range(0,len(phase)):
            phase[n] = (freq[n])*((time[n]-pulse_length/2)/1e3)*180
     
        phase_prog = np.zeros(int(nump))
        for n in range(0,len(phase)):
            phase_prog[n] = math.fmod(phase[n],360)  
       
    
    return time,amp,freq,phase_prog;