# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:47:58 2020

@author: philipp
"""


import numpy as np
import TD_NMR as td
import matplotlib.pyplot as plt
import lmfit as lmf
import math

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

def gauss(x,I0,T2):
    
    y = I0*np.exp(-0.5*(x/T2)**2)
    return y

def expo(x,I0,T2,n):
    y = I0*np.exp(-(x/T2)**n)
    return y


def abragam(x,I0,a,b):
    y = I0*np.exp(-0.5*(a*x)**2)*(np.sin(b*x))/(b*x)
    return y


def residual(params,x, data, noise):
    
    A1 = params['amp_1']
    A2 = params['a_1']
    # A3 = params['b_1']
    
    # B1 = params['amp_2']
    # B2 = params['T2_g']
    
    C1 = params['amp_3']
    C2 = params['T2_e']
    
    # D1 = params['noise']
    
    model = gauss(x,A1,A2)+expo(x,C1,C2,1)
    # +expo(x,B1,B2,1)
    # +gauss(x,B1,B2)+D1
    
    return np.sqrt((model-data)**2/len(x))

def final_model(params,x, data, noise):
    
    A1 = params['amp_1']
    A2 = params['a_1']
    # A3 = params['b_1']
    
    #B1 = params['amp_2']
   # B2 = params['T2_g']
    
    C1 = params['amp_3']
    C2 = params['T2_e']
    
    # D1 = params['noise']
    
    model = gauss(x,A1,A2)+expo(x,C1,C2,1)
    #+expo(x,B1,B2,1)
    # +gauss(x,B1,B2)+D1
    
    return model


#################################
## Test of MSE FID below
##
## full set of processing and fitting using LMFIT

#############

params = lmf.Parameters()

params.add('amp_1', value=50000)
params.add('a_1',value = 1e-6)
# params.add('b_1', value = 1, vary = False)

# params.add('amp_2', value = 5000)
# params.add('T2_g', value = 0.5e-3)

params.add('amp_3', value = 70000)
params.add('T2_e', value =1.4e-4)

# params.add('noise', value = 45)




test = td.TD_NMR_data('C:/Users/philipp/Box/NMR-data/500/200227-NIKE-MSE-test',221)
# correct for the digital filter oscillations and group delay points
test.correct_filter()
# shifts the time axis by subtracting the value
test.timeshift(10e-6)
test.truncate(2*1024)
# output values of the FID
# x2,y2 = test.return_fid()
x2,y2 =test.calc_mag()

print(x2)


model = residual(params,x2,y2,1)
out = lmf.minimize(residual, params, args=(x2,y2,1))
final_params = out.params
final = final_model(final_params,x2,y2,1)

print(final_params)
print(final)
fig = plt.figure()
# plt.scatter(x2,y2.imag)
# plt.scatter(x2,y2.real)
plt.scatter(x2,y2)
plt.plot(x2,final,color='red')
ax =plt.gca()
plt.xlim(-50e-6,2.5e-3)
plt.yscale('log')
plt.ylim(1e0,2e5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)
ax.set_xlabel(r'FID / s', fontdict=font)

plt.show()
# fig.savefig("NIKE-TDNMR.png",dpi=300, format='png', bbox_inches='tight', pad_inches=0.1)
# plt.close() 
