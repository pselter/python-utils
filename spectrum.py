# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 18:33:34 2019

@author: Selter
"""

from numba import njit, jit, complex128, float64
import numpy as np
import math
from scipy.linalg import expm,sinm,cosm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


##########################################################################################
##########################################################################################
##########################################################################################
### define rotation operations
### these give you the rotation matrix required for an operation

def adjunc(matrix):
    new_matrix = np.zeros(matrix.shape)
    new_matrix[:,0] = matrix[0,:]
    new_matrix[:,1] = matrix[1,:]
    new_matrix[:,2] = matrix[2,:]
    return new_matrix;

def x_rot(angle):
    theta = np.deg2rad(angle)
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta), -math.sin(theta) ],
                    [0,         math.sin(theta), math.cos(theta)  ]
                    ])
    #new_vec = R_x @ vec 
    return R_x;

def y_rot(angle):
    theta = np.deg2rad(angle)
    R_y = np.array([[math.cos(theta),    0,      math.sin(theta)  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta),   0,      math.cos(theta)  ]
                    ])
#    new_vec = R_y @ vec 
    return R_y;

def z_rot(angle):
    theta = np.deg2rad(angle)
    R_z = np.array([[math.cos(theta),    -math.sin(theta),    0],
                    [math.sin(theta),    math.cos(theta),     0],
                    [0,                     0,                      1]
                    ])
#    new_vec = R_z @ vec 
    return R_z;

def a_rot(axis,angle):
    u = axis
    norm = np.sqrt(u[0]**2+u[1]**2+u[2]**2)
    u = u/norm
    theta = np.deg2rad(angle)
    R_a = np.array([[math.cos(theta) + u[0]**2 * (1-math.cos(theta)), 
             u[0] * u[1] * (1-math.cos(theta)) - u[2] * math.sin(theta), 
             u[0] * u[2] * (1 - math.cos(theta)) + u[1] * math.sin(theta)],
            [u[0] * u[1] * (1-math.cos(theta)) + u[2] * math.sin(theta),
             math.cos(theta) + u[1]**2 * (1-math.cos(theta)),
             u[1] * u[2] * (1 - math.cos(theta)) - u[0] * math.sin(theta)],
            [u[0] * u[2] * (1-math.cos(theta)) - u[1] * math.sin(theta),
             u[1] * u[2] * (1-math.cos(theta)) + u[0] * math.sin(theta),
             math.cos(theta) + u[2]**2 * (1-math.cos(theta))]])
#    new_vec = R_a @ vec 
    return R_a;

##########################################################################################

##########################################################################################
# define how an ellipsoid is plotted    
# usefull for plotting the actual ellipsoid
        
def plot_CSA_ellipsoid(ax, tensor,  center,color='black',rstride=1, cstride=1,alpha=0.2,zorder=1000):
    """Plot the 3-d Ellipsoid ell on the Axes3D ax."""

    # points on unit sphere
    u = np.linspace(0.0, 2.0 * np.pi, 20)
    v = np.linspace(0.0, np.pi, 20)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    # transform points to ellipsoid
    for i in range(len(x)):
        for j in range(len(x)):
            x[i,j], y[i,j], z[i,j] = center+ np.dot(tensor,[x[i,j],y[i,j],z[i,j]])

    ax.plot_wireframe(x, y, z,  rstride=rstride, cstride=cstride, color=color, alpha=alpha,zorder=zorder)
                      
                      
def plot_CSA_axes(ax, tensor, center, color='red',alpha=0.8,zorder=0,normalize=True):
                     
    q0 = tensor[:,0]
    q1 = tensor[:,1]
    q2 = tensor[:,2]
    p0 = tensor[:,0]*-1.
    p1 = tensor[:,1]*-1.
    p2 = tensor[:,2]*-1.

    X, Y, Z = zip(center,center,center) 
    U, V, W = zip(q0,q1,q2)
    U2, V2, W2 = zip(p0,p1,p2)  
    
    ax.quiver(X,Y,Z,U,V,W,color=color,alpha=alpha,normalize=normalize,zorder=zorder)
##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################  
##########################################################################################
##########################################################################################
# Program to generate powder averaging according to the spherical ZCW method shown here:
# https://doi.org/10.1006/jmre.1996.1087
# Saving to crystal files not working correctly at this point
# you need to add the number of crystallites as the first line
# had no time to figure this out, maybe later
# REPULSION might follow at a later stage
 
def fibo(n): 
    """calculate the nth Fibonacci number in recursive way"""
    if n<0: 
        print("Incorrect input") 
    # First Fibonacci number is 0 
    elif n==0: 
        return 0
    # Second Fibonacci number is 1 
    elif n==1: 
        return 1
    else: 
        return fibo(n-1)+fibo(n-2) 

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def frac(n):
    trunc = truncate(n)
    return n-trunc


###################################################################################################  
# Discontinued, probably not correct:
# see: JOURNAL OF MAGNETIC RESONANCE 132, 220–239 (1998)
# "Computation of Orientational Averages in Solid-state NMR by Gaussian Spherical Quadrature"

#def zcw(N,file=False,deg=True,verbose=True):
#    """generate a zcw crystal file, with the Nth fibonacci number as the number of crystallites"""
#    Ncryst = fibo(N)
#    g1 = 1
#    g2 = fibo(N-2)
#    
#    alpha = 0
#    beta = 0
#    weight = 0 
#    #create the empty list
#    crysts = np.zeros((Ncryst,3))
#    if verbose == True:
#        print('creating crystal averaging values')
#        print('zcw scheme, number of crystallites: '+str(Ncryst))
#    if deg == False:
#        for a in range(1,Ncryst+1):
#            alpha = np.pi*2*frac(((a-1)/Ncryst)*g1)
#            beta = np.arccos(1-2*frac(((a-1)/Ncryst)*g2))
#            weight = np.sin(beta)
#            crysts[a-1,:] = alpha, beta, weight
#    else:
#        for a in range(1,Ncryst+1):
#            alpha = np.pi*2*frac(((a-1)/Ncryst)*g1)
#            beta = np.arccos(1-2*frac(((a-1)/Ncryst)*g2))
#            weight = np.sin(beta)
#            crysts[a-1,:] = np.rad2deg(alpha), np.rad2deg(beta), weight
#        
#    if file == False:
#        return crysts
#    else:
#       np.savetxt('zcw'+str(Ncryst)+'.cryst',crysts)
###################################################################################################
    
def STEP(N_beta,ab_ratio=4,file=False,deg=True,verbose=True):
    """generate a STEP crystal file, with orientations spanning a hemisphere"""
    
    Nb = N_beta
    Na = Nb*ab_ratio
    Ncryst = Nb*Na
    
    alphas = np.zeros(Na)
    betas = np.zeros(Nb)
    weight = 0
    stepweight = 0
    #create the empty list
    crysts = np.zeros((Ncryst,3))
    if verbose == True:
        print('creating crystal averaging values')
        print('STEP scheme, number of crystallites: '+str(Ncryst))
    
    for j in range(0,Nb):
        betas[j]=np.pi/(4*Nb)*(2*j+1)
        stepweight = stepweight + np.sin(betas[j])
        
    for i in range(0,Na):
        alphas[j]=(np.pi*2)/Na*i
        
            
    Nstep = 1/(Na*stepweight)
    
    k = 0
    if deg == False:
        for b in range(0,Nb):
            weight = Nstep*np.sin(betas[b])
            for a in range(0,Na): 
                crysts[k-1,:] = alphas[a], betas[b], weight
                k += 1
            
    else:
        for a in range(1,Ncryst+1):
           for b in range(0,Nb):
            weight = Nstep*np.sin(betas[b])
            for a in range(0,Na): 
                crysts[k-1,:] = np.rad2deg(alphas[a]), np.rad2deg(betas[b]), weight
                k += 1
        
    if file == False:
        return crysts
    else:
        np.savetxt('STEP'+str(Ncryst)+'.cryst',crysts)


def BETA(N_beta,file=False,deg=True,verbose=True):
    """generate a BETA crystal file, with only beta angles"""
    Ncryst = N_beta
    Nb = N_beta
    crysts = np.zeros((Ncryst,2))
    
    if deg == False:
        for j in range(0,Nb):
            beta = np.pi/(4*Nb)*(2*j+1)
            crysts[j,:]=beta,np.sin(beta)
    else:
        for j in range(0,Nb):
            beta = np.pi/(4*Nb)*(2*j+1)
            crysts[j,:]=np.rad2deg(beta),np.sin(beta)
    
    if verbose == True:
        print('creating Beta angle averaging values')
        print('number of beta angles: '+str(Ncryst))

    if file == False:
        return crysts
    else:
        np.savetxt('beta'+str(Ncryst)+'.cryst',crysts)
        
     

def ZCW(M,ab_ratio=4,file=False,deg=True,verbose=True):
    """generate a zcw crystal file, with orientations spanning a hemisphere"""
    
    M = M+6
        
    Ncryst = fibo(M+2)
    g_m = fibo(M)
    
    alpha = 0
    beta = 0
    weight = 0 
    #create the empty list
    crysts = np.zeros((Ncryst,3))
    if verbose == True:
        print('creating crystal averaging values')
        print('zcw scheme, number of crystallites: '+str(Ncryst))
    if deg == False:
        for j in range(0,Ncryst):
            alpha = (np.pi*2)*math.fmod((j*g_m/Ncryst),1)
            beta = np.arccos(2*math.fmod(j/Ncryst,1)-1)
            weight = 1/Ncryst
            crysts[j-1,:] = alpha, beta, weight
    else:
        for j in range(0,Ncryst):
            alpha = (np.pi*2)*math.fmod((j*g_m/Ncryst),1)
            beta = np.arccos(2*math.fmod(j/Ncryst,1)-1)
            weight = 1/Ncryst
            crysts[j-1,:] = np.rad2deg(alpha), np.rad2deg(beta), weight
   
    if file == False:
        return crysts
    else:
        np.savetxt('ZCW'+str(Ncryst)+'.cryst',crysts)
        
##########################################################################################
# This calculates the frequency in the lab frame
# omega0 is the carrier frequency
# tensor is a tensor object
# alpha, beta, and gamma are the three euler angles     
        
def lab_cs(omega0,tensor,alpha=0,beta=0,gamma=0):
    
    """this calculates the resonance frequency
    using an arbitrary csa tensor and its relative orientation
    
    - requires the tensor property of a tensor class object
    - omega_0 is the larmor frequency of the nucleus in question in Hz
    - alpha, beta, gamma are the three euler angles    
    - the elements of the tensor are assumed to be in ppm
    """
    
    
    # define column vector
    col_vec = np.array([[0],[0],[1]])
    field_vector = np.array([0,0,1])
    
    # rotate the csa tensor from Principle axes frame to the lab frame
    alpha_rot = z_rot(alpha)
    adjunc_alpha_rot  = adjunc(alpha_rot)
    beta_rot = y_rot(beta)
    adjunc_beta_rot  = adjunc(beta_rot)
    gamma_rot = z_rot(gamma)
    adjunc_gamma_rot  = adjunc(gamma_rot)

    rot_tensor= adjunc_gamma_rot @ adjunc_beta_rot @ adjunc_alpha_rot @ tensor @ alpha_rot @ beta_rot @ gamma_rot
    
    vector_prod = field_vector @ rot_tensor @ col_vec
    omega_cs = (-1.0*omega0*(vector_prod*1e-6))

    return omega_cs
##########################################################################################


def single_cryst_signal(larmor_freq,tensor,alpha,beta,timeaxis):
    
    omega1=lab_cs(larmor_freq,tensor,alpha,beta)
    signal = np.exp(-1j*omega1*timeaxis)
    
    return signal

def single_cryst_echo(larmor_freq,tensor,alpha,beta,timeaxis,tau=1e-4):
    
    omega1=lab_cs(larmor_freq,tensor,alpha,beta)
    
    echo_signal = np.exp(1j*omega1*tau)
    signal = np.exp(-1j*omega1*(timeaxis-tau))*echo_signal
    
    return signal


def powder_signal(larmor_freq,tensor,crystal_list,timeaxis,T2=1e-3):
    
    n_cryst = len(crystal_list[:,0])
    
    
    avg_signal = single_cryst_signal(larmor_freq,tensor,crystal_list[0,0],crystal_list[0,1],timeaxis)*np.exp(-1.0*timeaxis/T2)
    
    for n in range(1,n_cryst):
        avg_signal += single_cryst_signal(larmor_freq,tensor,crystal_list[n,0],crystal_list[n,1],timeaxis)*np.exp(-1.0*timeaxis/T2)
    
   
    return avg_signal



def two_site_echo(larmor_freq,tensor_A,tensor_B,alpha,beta,timeaxis,tau=1e-4,T_exchange=1e-4):
    
    omega1=lab_cs(larmor_freq,tensor_A,alpha,beta)
    omega2=lab_cs(larmor_freq,tensor_B,alpha,beta)
    
    
    echo_signal1 = np.exp(1j*omega1*tau)*np.exp(-1.0*tau/T_exchange)
   
    echo_signal2 = np.exp(1j*omega2*tau)*(1-np.exp(-1.0*tau/T_exchange))

    
#    signal = (np.exp(-1j*omega1*(timeaxis-tau))*echo_signal1*np.exp(-1.0*(timeaxis-tau)/T_exchange))+(np.exp(-1j*omega2*(timeaxis-tau))*echo_signal2*(1-np.exp(-1.0*(timeaxis-tau)/T_exchange)))+(np.exp(-1j*omega1*(timeaxis-tau))*echo_signal2*np.exp(-1.0*(timeaxis-tau)/T_exchange))+(np.exp(-1j*omega2*(timeaxis-tau))*echo_signal1*(1-np.exp(-1.0*(timeaxis-tau)/T_exchange)))
#    
    return signal
    
    
    
    
def powder_2site_echo(larmor_freq,tensor_A,tensor_B,crystal_list,timeaxis,tau=1e-4,T2=1e-3,T_exchange=1e-4):
    
    n_cryst = len(crystal_list[:,0])
    
    
    avg_signal = two_site_echo(larmor_freq,tensor_A,tensor_B,crystal_list[0,0],crystal_list[0,1],timeaxis,tau,T_exchange)*np.exp(-1.0*timeaxis/T2)
    
    for n in range(1,n_cryst):
        avg_signal += two_site_echo(larmor_freq,tensor_A,tensor_B,crystal_list[n,0],crystal_list[n,1],timeaxis,tau,T_exchange)*np.exp(-1.0*timeaxis/T2)
    
   
    return avg_signal


##########################################################################################







##########################################################################################
##########################################################################################
##########################################################################################        





##########################################################################################
##########################################################################################
##########################################################################################
#
# Definition of the tensor class
# usefull for handling 3x3 tensors/matrices
# WIP, use with caution, duh
# mainly used to convert stuff and to plot the damn thing

class generic_tensor(object):
    
    """generic tensor object
    follows the Haeberlen convention
    
    zz,yy,xx are the diagonal components
    can give you the following values:
            
    isotropic d_iso
    anisotropy          R_zz - (R_xx + R_yy)
    reduced anisotropy  R_zz - R_iso
    eta                 (R_yy - R_xx) / (R_zz - R_iso)
    
    beware the confusion with anisotropy and the reduced anisotropy!!
    """  
    
    def __init__(self,iso,delta,eta):
        
        self.iso = iso
        self.red_aniso = delta
        self.eta = eta
        
        
        if delta > 0:
            self.d11 = self.iso+self.red_aniso
            self.d22 = self.iso-self.red_aniso*(1-self.eta)/2
            self.d33 = self.iso-self.red_aniso*(1+self.eta)/2
            
            self.dzz = self.d11
            self.dyy = self.d22
            self.dxx = self.d33
        else:
            self.d33 = self.iso+self.red_aniso
            self.d22 = self.iso-self.red_aniso*(1-self.eta)/2
            self.d11 = self.iso-self.red_aniso*(1+self.eta)/2

            self.dzz = self.d33
            self.dyy = self.d22
            self.dxx = self.d11
            
        self.tensor =  np.array([[self.dxx,0,0],[0,self.dyy,0],[0,0,self.dzz]])
        self.Rtensor =  np.array([[self.dxx-iso,0,0],[0,self.dyy-iso,0],[0,0,self.dzz-iso]])
        self.aniso = 3*delta/2
#        
    def iso(self):
        zolo = self.iso
        return zolo
    
    def aniso(self):
        zdelta = self.aniso
        return zdelta
    
    def red_aniso(self):
        zred_aniso = self.red_aniso
        return zred_aniso
    
    def delta(self):
        zred_aniso = self.red_aniso
        return zred_aniso
    
    def eta(self):
        zeta = self.eta
        return zeta

    def tensor(self):
        ztensor = self.tensor
        return ztensor
    
    def plotmeR(self):
        fig = plt.figure(figsize=(4., 4.),facecolor=None, frameon=False)
        ax = fig.add_subplot(111, projection='3d')
        ax = plt.gca(projection='3d')
#       ax._axis3don = False
        
        plot_CSA_axes(ax,self.Rtensor,(0,0,0),normalize=False,color='blue')
        plot_CSA_ellipsoid(ax,self.Rtensor,(0,0,0),color='blue')
        tensor = self.Rtensor
        ax.set_xlim(-1*tensor[2,2],tensor[2,2])
        ax.set_ylim(-1*tensor[2,2],tensor[2,2])
        ax.set_zlim(-1*tensor[2,2],tensor[2,2])
        ax.set_aspect('equal')
        fig.tight_layout()
        ax.view_init(20, 300)
        plt.show()

        plt.close()
        
    def plotme(self):
        fig = plt.figure(figsize=(4., 4.),facecolor=None, frameon=False)
        ax = fig.add_subplot(111, projection='3d')
        ax = plt.gca(projection='3d')
#       ax._axis3don = False
        tensor = self.tensor
        plot_CSA_axes(ax,tensor,(0,0,0),normalize=False,color='blue')
        plot_CSA_ellipsoid(ax,tensor,(0,0,0),color='blue')
        
        ax.set_xlim(-1*tensor[1,1],tensor[1,1])
        ax.set_ylim(-1*tensor[0,0],tensor[0,0])
        ax.set_zlim(-1*tensor[2,2],tensor[2,2])
        ax.set_aspect('equal')
        fig.tight_layout()
        ax.view_init(20, 300)
        plt.show()

        plt.close()
        
### END OF DEFINITIONS
##########################################################################################
##########################################################################################
##########################################################################################        






























##########################################################################################
# Define functions to generate a powder pattern
# check carefully whether the definitions of delta and eta are correct!!!   
# this does not work in the tensor object for some f***** up reason   
#
# Using numba this is actually pretty fast    

#############################################
# This is the function that generates the frequency of a given crystallite
@njit
def omega(theta,phi,delta,eta):
    omega = (1/2)*delta*(3*(np.cos(theta)**2)-1-eta*np.sin(theta)*np.cos(2*phi))
    return omega


#############################################
# This is the function that generates the powder pattern of a static sample
# cristallites are staying in their relative positions, so no time dependence here
# calculates the intensities in the frequency domain
@njit
def powder_spectrum(axis,iso=0.0,delta=5.0,eta=0.5,res=1.0,width=1.0):

    deg_thetas = np.arange(0,180,res)
    thetas = np.deg2rad(deg_thetas)
    deg_phis = np.arange(0,360,res)
    phis = np.deg2rad(deg_phis)
    intensity=np.zeros(len(axis))
    
    for k in thetas:
        for l in phis:
#            omega = iso-(1/2)*delta*(3*(np.cos(k)**2)-1-eta*np.sin(k)*np.cos(2*l))
            total_omega = iso-omega(k,l,delta,eta)
            curint = np.exp(-np.power(axis - total_omega, 2.0) / (2 * np.power(width, 2.0)))*np.sin(k)
            intensity += curint
    return intensity 
##########################################################################################

