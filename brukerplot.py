import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt

class bruker1d(object):
		
	def __init__(self, path, expno=1, procno=1):
		
		self.dic, self.data = ng.bruker.read_pdata(path+'/'+str(expno)+'/pdata/'+str(procno))
		self.pdata = ng.bruker.read_procs_file(path+'/'+str(expno)+'/pdata/'+str(procno))
		self.udic = ng.bruker.guess_udic(self.dic,self.data)
		self.uc = ng.fileiobase.uc_from_udic(self.udic, 0)
		#print('bruker object created')
		
	def plot1d(self):
		self.xaxis = self.uc.ppm_scale()
		self.yaxis = self.data
	#	self.xaxis = np.linspace(start,stop, num=len(self.data[self.uc(str(start)+' ppm'): self.uc(str(stop)+' ppm')]))
	#	self.yaxis = self.data[self.uc(str(start)+' ppm'): self.uc(str(stop)+' ppm')]
		return self.xaxis, self.yaxis
		
		
		
class bruker2d(object):
	
	
    def __init__(self, path, expno=1, procno=1):
	
        self.dic, self.data = ng.bruker.read_pdata(path+'/'+str(expno)+'/pdata/'+str(procno))
        self.pdata = ng.bruker.read_procs_file(path+'/'+str(expno)+'/pdata/'+str(procno))
        self.udic = ng.bruker.guess_udic(self.dic,self.data)
        self.uc0 = ng.fileiobase.uc_from_udic(self.udic, 0)
        self.uc1 = ng.fileiobase.uc_from_udic(self.udic, 1)
        #print('bruker object created')
		
    def plot2d(self):
        x_values = self.uc1.ppm_scale()
        y_values = self.uc0.ppm_scale()
        self.x, self.y = np.meshgrid(x_values,y_values)
        return self.x, self.y, self.data
        
    def plot_rocsa(self):
        x_values = self.uc1.ppm_scale()
        y_values = self.uc0.hz_scale()
        self.x, self.y = np.meshgrid(x_values,y_values)
        return self.x, self.y, self.data
		
    def contour_level(self,zoom=1.0, mult=1.8, n_pos_levels=8, n_neg_levels=8):
	
		#contours, do not mess with this

        cl = (self.data.std() * zoom * mult ** np.arange(n_pos_levels))
        cl2 = (self.data.std() * (-1.0*zoom) * mult ** np.arange(n_neg_levels))
        cl3 = cl2[::-1]
        self.c_levels = [*cl3, *cl]
        return self.c_levels
		# contours end

    def proj_1(self, start, stop):	
		
        self.proj_1_scale = self.uc1.ppm_scale()
        self.proj_1_data = np.arange(0,len(self.data[0,0:]),1)
				
        for k in self.proj_1_data:
            test = sum(x for x in self.data[self.uc0(str(start)+' ppm'):self.uc0(str(stop)+' ppm'),k] if x > 0)
            self.proj_1_data[k]=test
		
        return self.proj_1_scale, self.proj_1_data
    
    def slice_0_hz(self,line):
        self.slice_0_scale = self.uc0.hz_scale()
        self.slice_0_data = self.data[0:,line]
				  		
        return self.slice_0_scale, self.slice_0_data 
	
    def slice_0(self,line):
        self.slice_0_scale = self.uc0.ppm_scale()
        self.slice_0_data = self.data[0:,line]
				  		
        return self.slice_0_scale, self.slice_0_data 
		
		
    def proj_0(self, start, stop):
		
        self.proj_0_scale = self.uc0.ppm_scale()
        self.proj_0_data = np.arange(0,len(self.data[0:,0]),1)
				
        for k in self.proj_0_data:
            test = sum(x for x in self.data[k,self.uc1(str(start)+' ppm'):self.uc1(str(stop)+' ppm')] if x > 0)
            self.proj_0_data[k]=test
		
        return self.proj_0_scale, self.proj_0_data
    
    def proj_0_hz(self, start, stop):
		
        self.proj_0_scale = self.uc0.hz_scale()
        self.proj_0_data = np.arange(0,len(self.data[0:,0]),1)
				
        for k in self.proj_0_data:
            test = sum(x for x in self.data[k,self.uc1(str(start)+' ppm'):self.uc1(str(stop)+' ppm')] if x > 0)
            self.proj_0_data[k]=test
		
        return self.proj_0_scale, self.proj_0_data
		
		