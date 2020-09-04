#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 16:50:59 2018

@author: akhansen
"""

import numpy as np
import math


#################################
### define rotation operations
def x_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta), -math.sin(theta) ],
                    [0,         math.sin(theta), math.cos(theta)  ]
                    ])
    new_vec = np.dot(R_x,vec)
    return new_vec;

def y_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_y = np.array([[math.cos(theta),    0,      math.sin(theta)  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta),   0,      math.cos(theta)  ]
                    ])
    new_vec = np.dot(R_y,vec)
    return new_vec;

def z_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_z = np.array([[math.cos(theta),    -math.sin(theta),    0],
                    [math.sin(theta),    math.cos(theta),     0],
                    [0,                     0,                      1]
                    ])
    new_vec = np.dot(R_z,vec)
    return new_vec;

def a_rot(vec,axis,angle):
    u = axis
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
    new_vec = np.dot(R_a,vec)
    return new_vec;
#################################
#################################

####in which frame of reference do you want this calculation:
case = 'rotframe'
#case = 'freqmodframe'


### strength frequency sweep 
f_start =  0.2
f_end   = -0.2

### Relative strength of the B1 field
b1_nutation = 0.001

### Amplitude factor for shape
N = 80

## b1_nutation = 0.025 excitation
## b1_nutation = 0.1   inversion 
#################################

#time axis
t_start = 0
t_stop = 300
t_points = 5000

### magnetization vector: x,y,z,time
start_vec = np.array([0.0,0.0,1,0])
b1_vec = np.array([0,-1,0])
ofs_vec = np.array([0,0,1])

time = np.linspace(t_start,t_stop,t_points)
t_step = (t_stop-t_start)/t_points

f_diff = f_end-f_start
f_step = f_diff/len(time)
freq = np.zeros(len(time))
for n in range(0,len(freq)):
    freq[n]=f_start+n*f_step

#################################
#### calculate the frequency ramp    
#amp = b1_nutation*1.0+time*0.0
amp = b1_nutation*(1-np.abs(np.cos((np.pi*time)/(t_stop-t_start))**N))
#################################

print('freq rate is: '+str(f_step/t_step))
print('freq_r/b1 ratio is '+str((f_step/t_step)/b1_nutation))
print('b1 ratio/freq_r is '+str(b1_nutation/(f_step/t_step)))

nut_angle = amp*t_step*360
mag_vec = np.zeros((len(time),4))
b1_vect = np.zeros((len(time),3))
beff_vec = np.zeros((len(time),4))


beff_vec[:,3] = time
beff_vec[:,2] = freq*-1.0
beff_vec[:,:1] = b1_vec[:1]
beff_vec[:,0] = beff_vec[:,0]*amp*-1.0
beff_vec[:,1] = beff_vec[:,1]*amp*-1.0

mag_vec[0] = start_vec
b1_vect[0] = b1_vec
offset_angle = freq*t_step*360
arcos_vec = np.zeros((len(time),4))
arcos_vec[:,3] = time
betrag = np.zeros((len(time),1))
theta = np.zeros((len(time),1))
ratio = np.zeros((len(time),1))
ratio[0] = 1

if case == 'rotframe':
    for n in range(1,len(time)):
        mag_vec[n,:3] = a_rot(mag_vec[n-1,:3],b1_vect[n-1],nut_angle[n])
        mag_vec[n,3] = time[n]
        b1_vect[n] = z_rot(b1_vect[n-1],offset_angle[n])
        ratio[n]=offset_angle[n]/nut_angle[n]
        
        
if case == 'freqmodframe':
    for n in range(1,len(time)):
        mag_vec[n,:3] = a_rot(mag_vec[n-1,:3],b1_vect[n-1],nut_angle[n])
        mag_vec[n,:3] = z_rot(mag_vec[n,:3],offset_angle[n])
        mag_vec[n,3] = time[n]
        b1_vect[n] = b1_vect[n-1]
        arcos_vec[n,2] = math.acos(mag_vec[n,2])
        betrag[n] = math.sqrt(freq[n]**2+amp[n]**2)
        theta[n] = math.acos(freq[n]/betrag[n])
        ratio[n]=nut_angle[n]-offset_angle[n]

import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


xx = mag_vec[:,0]
yy = mag_vec[:,1]
zz = mag_vec[:,2]
tt = mag_vec[:,3]
xxx =  arcos_vec[:,0]
yyy =  arcos_vec[:,1]
zzz =  arcos_vec[:,2]
ttt =  arcos_vec[:,3]

#Set colours and render
fig = plt.figure(dpi=150,figsize=(4,4))
ax = plt.axes(projection='3d')





plotter = ax.scatter3D(xx,yy,zz,c=tt,zdir='z',s=10,cmap='plasma',zorder=10)
#ax.scatter(xxx,yyy,zzz,c=ttt,s=10,cmap='seismic')


# draw sphere
ua, va = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
xa = np.cos(ua)*np.sin(va)
ya = np.sin(ua)*np.sin(va)
za = np.cos(va)
ax.plot_wireframe(xa, ya, za, color="k",linewidth=0.5,alpha=0.5)
# ax.set_axisbelow(True)


ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
# ax.set_aspect("equal")
plt.tight_layout()
plt.show()


fig2 = plt.figure()
#plt.plot(time,np.absolute(offset_angle))
#plt.plot(time,nut_angle)
#plt.plot(time,-1*nut_angle)
plt.plot(time,mag_vec[:,2],label='z mag')
plt.plot(time,mag_vec[:,1],label='y mag')
plt.plot(time,mag_vec[:,0],label='x mag')
plt.legend()
#plt.plot(time,arcos_vec[:,2])
#plt.plot(time,theta)
#plt.plot(time,ratio)
plt.ylim(-1,1)
plt.show()


#
#fig3 = plt.figure()
#plt.plot(time,mag_vec[:,0])
#plt.plot(time,mag_vec[:,1])
#plt.plot(time,mag_vec[:,2])
#
##plt.show()

