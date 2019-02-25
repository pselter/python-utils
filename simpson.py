# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 15:07:00 2019

@author: Selter
"""

###############
#INFO:
#
# simple little code designed to provide easy SIMPSON access in python
# adapt the location of your simpson executable if necessary
#
# current aim is to just use it for starting simulations
# inputfile handling is done by different program, butmight be included in the future
##############


import numpy as np
import subprocess

### SIMPSON program path, only change if you have this installed somewhere else
simpson = 'C:\Program Files (x86)\SIMPSON\simpson.exe'


    
def py2spe(real,imag,file,np,sw,typ='SPE',ref=0.0):
    header = 'SIMP'+'\n'+'NP='+str(np)+'\n'+'SW='+str(sw)+'\n'+'REF='+str(ref)+'\n'+'TYPE='+str(typ)+'\n'+'DATA'+'\n'
    
    spefile = open(file, 'w', newline='\n')
    spefile.write(header)
    for n in range(len(real)):
        spefile.write(str(real[n])+' '+str(imag[n])+'\n')        
    spefile.write('END '+'\n')
    spefile.close()
     


class sim(object):
    
    def __init__(self,inputfile):
        self.exec = simpson
        self.infile=inputfile
       
    def runsim(self,output=True):
        self.prog = subprocess.run([self.exec,self.infile],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
        
        if output==True:
            return self.prog.returncode, self.prog.stdout, self.prog.stderr
            
        return 
        
#### TO DO
###  inputfile generation
        # inputfile manipulation
        # changing calculation parameters????