# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 22:37:10 2018

@author: philipp
"""

from numba import jit
import numpy as np
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



#################################
### define rotation operations

def adjunc(matrix):
    new_matrix = np.zeros(matrix.shape)
    new_matrix[:,0] = matrix[0,:]
    new_matrix[:,1] = matrix[1,:]
    new_matrix[:,2] = matrix[2,:]
    return new_matrix;

def x_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta), -math.sin(theta) ],
                    [0,         math.sin(theta), math.cos(theta)  ]
                    ])
    #new_vec = R_x @ vec 
    return R_x;

def y_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_y = np.array([[math.cos(theta),    0,      math.sin(theta)  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta),   0,      math.cos(theta)  ]
                    ])
#    new_vec = R_y @ vec 
    return R_y;

def z_rot(vec,angle):
    theta = np.deg2rad(angle)
    R_z = np.array([[math.cos(theta),    -math.sin(theta),    0],
                    [math.sin(theta),    math.cos(theta),     0],
                    [0,                     0,                      1]
                    ])
#    new_vec = R_z @ vec 
    return R_z;

def a_rot(vec,axis,angle):
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

#################################
    

#################################
# define how an ellipsoid is plotted    
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
    
    print('yolo this is')
    print(W)
#    for n in range(3):
#        ax.plot((U[n],center[n]),(V[n],center[n]),(W[n],center[n]),color='black',alpha=0.2)
#        print((U[n],center[n]))
#        ax.plot((U2[n],center[n]),(V2[n],center[n]),(W2[n],center[n]),color='black',alpha=0.2)
#     
    ax.quiver(X,Y,Z,U,V,W,color=color,alpha=alpha,normalize=normalize,zorder=zorder)
  
#################################



#################################
# define the test tensor
center0= (0,0,0)
center2= (-1.5,0,0)    
center3= (0,0,0) 
center4= (1.5,0,0) 
center5= (3,0,0) 
sigma = np.array([['xx','xy','xz'],['yx','yy','yz'],['zx','zy','zz']])
print(sigma)

sigma = np.zeros((3,3))
unity = np.zeros((3,3))
#original tensor is defined here
sigma[0,0]=0.2
sigma[1,1]=0.3
sigma[2,2]=.75


#unity tensor is defined here
unity[0,0]=1
unity[1,1]=1
unity[2,2]=1

#rotation matrix 1
rot_mat = a_rot(unity,([1,0,0]),60)
rot_mat2 = a_rot(unity,([1,0,0]),-30)
##################################
## FIGURE #1
rot_sigma = rot_mat @  sigma

transformed = rot_mat @ sigma @ rot_mat.transpose()
transformed_eigenvalues, transformed_eigenvectors = np.linalg.eig(transformed)

print("this is sigma")
print(sigma)

print("this is rotated sigma")
print(rot_sigma)

print("this is transformed")
print(transformed)


average = (sigma*0.5 + transformed*0.5) 
print("this is the average")

average_eigenvalues,average_eigenvectors = np.linalg.eig(average)
print(average)
yolo = average_eigenvectors
print(yolo)
#average_pas = 

#
#print("this is the average in PAS?")
#print(diagonal)

new_vec = np.zeros((3,3))
new_vec[:,0] = average_eigenvectors[:,2]
new_vec[:,1] = average_eigenvectors[:,1]*-1 
new_vec[:,2] = average_eigenvectors[:,0]*-1

print(new_vec)

average_plot = average @ new_vec
print(average_plot)
###############################


#################################
#### FIGURE #2
#trans_sigma = rot_mat @ sigma @ rot_mat.transpose()
#rot_sigma = rot_mat @ sigma
#print(trans_sigma)
#print(rot_sigma)
#eigen = rot_mat @ unity
#
#


###############################
# plot for standard tensor
###############
fig = plt.figure(figsize=(10., 10.),facecolor=None, frameon=False)
ax = fig.add_subplot(111, projection='3d')
ax = plt.gca(projection='3d')
ax._axis3don = False

plot_CSA_axes(ax,sigma,center0,normalize=False,color='blue')
plot_CSA_ellipsoid(ax,sigma,center0,color='blue')
print(sigma)

ax.set_xlim(-1, 1)
ax.set_ylim(1, -1)
ax.set_zlim(-1,1)
ax.set_aspect('equal')
fig.tight_layout()
ax.view_init(20, 300)
#for angle in range(0, 360):
#    ax.view_init(30, angle)
#    plt.draw()
#    plt.pause(.001)
plt.show()
#fig.savefig("Tensor_test4.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
plt.close()
#################################


##################################
##plot the rotated tensor
#################################
#fig = plt.figure(figsize=(10., 10.),facecolor=None, frameon=False)
#ax = fig.add_subplot(111, projection='3d')
#ax = plt.gca(projection='3d')
#ax._axis3don = False
#
#
#plot_CSA_axes(ax,rot_sigma,center0,normalize=False,color='red')
#plot_CSA_ellipsoid(ax,rot_sigma,center0,color='red')
#print(rot_sigma)
#
#
#ax.set_xlim(-1, 1)
#ax.set_ylim(1, -1)
#ax.set_zlim(-1,1)
#ax.set_aspect('equal')
#fig.tight_layout()
#ax.view_init(20, 300)
##for angle in range(0, 360):
##    ax.view_init(30, angle)
##    plt.draw()
##    plt.pause(.001)
#plt.show()
#fig.savefig("Tensor_test2.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
#plt.close()
##
#
#

##################################
##plot the rotated tensor in the CAF
#################################
#fig = plt.figure(figsize=(10., 10.),facecolor=None, frameon=False)
#ax = fig.add_subplot(111, projection='3d')
#ax = plt.gca(projection='3d')
#ax._axis3don = False
#
#
#plot_CSA_axes(ax,transformed,center0,normalize=False,color='blue')
#plot_CSA_ellipsoid(ax,transformed,center0,color='blue')
#print(rot_sigma)
#
#
#ax.set_xlim(-1, 1)
#ax.set_ylim(1, -1)
#ax.set_zlim(-1,1)
#ax.set_aspect('equal')
#fig.tight_layout()
#ax.view_init(20, 300)
##for angle in range(0, 360):
##    ax.view_init(30, angle)
##    plt.draw()
##    plt.pause(.001)
#plt.show()
#fig.savefig("Tensor_test3.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
#plt.close()
##
##################################
##plot the average tensor in the CAF
#################################
#fig = plt.figure(figsize=(10., 10.),facecolor=None, frameon=False)
#ax = fig.add_subplot(111, projection='3d')
#ax = plt.gca(projection='3d')
#ax._axis3don = False
#
#
#plot_CSA_axes(ax,average,center0,normalize=False,color='blue')
#plot_CSA_ellipsoid(ax,average,center0,color='blue')
#print(average)
#
#
#ax.set_xlim(-1, 1)
#ax.set_ylim(1, -1)
#ax.set_zlim(-1,1)
#ax.set_aspect('equal')
#fig.tight_layout()
#ax.view_init(20, 300)
##for angle in range(0, 360):
##    ax.view_init(30, angle)
##    plt.draw()
##    plt.pause(.001)
#plt.show()
#fig.savefig("Tensor_test5.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
#plt.close()

##################################
##plot the average tensor in its PAF
#################################
#fig = plt.figure(figsize=(10., 10.),facecolor=None, frameon=False)
#ax = fig.add_subplot(111, projection='3d')
#ax = plt.gca(projection='3d')
#ax._axis3don = False
#
#
#eigen_norm = unity
#eigen_av = average_eigenvectors * average_eigenvalues
#eigen_rot = transformed_eigenvectors * transformed_eigenvalues
#
#
#plot_CSA_ellipsoid(ax,average_plot,center0,color='green')
#plot_CSA_axes(ax,eigen_av,center3,color='green',normalize=False,alpha=0.9)
#
#print(eigen_av)
#
#
#ax.set_xlim(-1, 1)
#ax.set_ylim(1, -1)
#ax.set_zlim(-1,1)
#ax.set_aspect('equal')
#fig.tight_layout()
#ax.view_init(20, 300)
##for angle in range(0, 360):
##    ax.view_init(30, angle)
##    plt.draw()
##    plt.pause(.001)
#plt.show()
#fig.savefig("Tensor_test6.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
#plt.close()

















