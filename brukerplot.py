import numpy as np
import nmrglue as ng
import matplotlib.pyplot as plt


class bruker1d(object):
		
    def __init__(self, path, expno=1, procno=1):
        self.dic, self.data = ng.bruker.read_pdata(path+'/'+str(expno)+'/pdata/'+str(procno))
        self.pdata = ng.bruker.read_procs_file(path+'/'+str(expno)+'/pdata/'+str(procno))
        self.udic = ng.bruker.guess_udic(self.dic,self.data)
        self.uc = ng.fileiobase.uc_from_udic(self.udic, 0)
           
    def plot1d(self):
        self.xaxis = self.uc.ppm_scale()
        self.yaxis = self.data
        return self.xaxis, self.yaxis


    def plot1d_norm(self):
        
        self.xaxis = self.uc.ppm_scale()
        self.yaxis = np.array(self.data)
        self.normaxis = self.yaxis/self.yaxis.max()
	#	self.xaxis = np.linspace(start,stop, num=len(self.data[self.uc(str(start)+' ppm'): self.uc(str(stop)+' ppm')]))
	#	self.yaxis = self.data[self.uc(str(start)+' ppm'): self.uc(str(stop)+' ppm')]
        return self.xaxis, self.normaxis
    
    
    def quickplot(self,xlim=(200,0),name='quickplot.png',savefig=True):
        
        self.font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }
       
        self.x, self.y  = self.plot1d_norm()
        
        self.fig = plt.figure(dpi=150, figsize=(6,4))
        self.ax = plt.gca()
        
        plt.plot(self.x, self.y,'-',color='k',linewidth=1)
        
        
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelsize=10)
        plt.tick_params(axis='y', which='both', labelleft=False, left=False)
        plt.minorticks_on()
        plt.xlim(xlim)
        self.ax.set_xlabel(r'$\delta_{iso}$'+' / ppm', fontdict=self.font)
        plt.tight_layout()
        plt.show()
        # save the figure
        if savefig == True:
            self.fig.savefig(name,dpi=300, format='png', bbox_inches='tight', pad_inches=0.1) 
        plt.close()
        return self.fig, self.ax
    
    def calc_zfratio(self):
        self.tdeff = self.dic["procs"]["TDeff"]
        self.si = self.dic["procs"]["SI"]
        self.zf_fac = self.si/self.tdeff
        self.sn_fac = np.sqrt(1/self.zf_fac)
        return self.zf_fac , self.sn_fac
        	
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
        
  
    def def_contours(self,zoom=1.0, mult=1.8, n_pos_levels=8, n_neg_levels=8):
	
		#contours, do not mess with this

        cl = (self.data.std() * zoom * mult ** np.arange(n_pos_levels))
        cl2 = (self.data.std() * (-1.0*zoom) * mult ** np.arange(n_neg_levels))
        cl3 = cl2[::-1]
        self.c_levels = [*cl3, *cl]
#        return self.c_levels
		# contours end
        
    def def_contours2(self,zoom=1.0, mult=1.8, n_pos_levels=8, n_neg_levels=8):
	
		#contours, do not mess with this
        # print(self.data.std())
        cl = ( zoom * mult ** np.arange(n_pos_levels))
        cl2 = ( (-1.0*zoom) * mult ** np.arange(n_neg_levels))
        cl3 = cl2[::-1]
        self.c_levels = [*cl3, *cl]
#        return self.c_levels
		# contours end
        
    def def_contours_percent(self,zoom=1.0, inc=10.0, n_pos_levels=8, n_neg_levels=8):
	
		#contours, do not mess with this

        cl = ( zoom + inc * np.arange(n_pos_levels))
        cl2 = ( (-1.0*zoom)  - inc * np.arange(n_neg_levels))
        cl3 = cl2[::-1]
        self.c_levels = [*cl3, *cl]
#        return self.c_levels
		# contours end

    def proj_1(self, start, stop):	
		
        self.proj_1_scale = self.uc1.ppm_scale()
        self.proj_1_data = np.arange(0,len(self.data[0,0:]),1)
				
        for k in self.proj_1_data:
            test = sum(x for x in self.data[self.uc0(str(start)+' ppm'):self.uc0(str(stop)+' ppm'),k] if x > 0)
            self.proj_1_data[k]=test
		
        return self.proj_1_scale, self.proj_1_data

    
    def proj_1_rocsa_norm(self, start, stop):	
		
        self.proj_1_scale = self.uc1.ppm_scale()
        self.proj_1_data = np.arange(0,len(self.data[0,0:]),1)
    
        self.ndata = self.data/self.data.max()*100
        
        for k in self.proj_1_data:
            test = sum(x for x in self.ndata[:,k] if x > 0)
            self.proj_1_data[k]=test
		
       
        return self.proj_1_scale, self.proj_1_data


    def proj_1_2(self, start, stop):	
        # if start < stop :
        #     store = start
        #     start = stop
        #     stop = store
            
        self.proj_1_scale = self.uc1.ppm_scale()
        self.proj_1_data = np.arange(0,len(self.data[0,0:]),1)
        print(self.uc0(str(start)+' ppm'))
        print(self.uc0(str(stop)+' ppm'))
		
        for k in self.proj_1_data:
            test = sum(x for x in self.data[self.uc0(str(start)+' ppm'):self.uc0(str(stop)+' ppm'),k] if x > 0)
            self.proj_1_data[k]=test
            
        self.normproj1 = self.proj_1_data
        self.nproj1 = np.array(self.normproj1)
        # print(self.nproj1)
        
        self.normproj12 = self.nproj1/self.nproj1.max()*100
        
        return self.proj_1_scale, self.normproj12

    
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


    def proj_0_2(self, start, stop):
		
        self.proj_0_scale = self.uc0.ppm_scale()
        self.proj_0_data = np.arange(0,len(self.data[0:,0]),1)
		
        for k in self.proj_0_data:
            test = sum(x for x in self.data[k,self.uc1(str(start)+' ppm'):self.uc1(str(stop)+' ppm')] if x > 0)
            self.proj_0_data[k]=test

        self.normproj0 = self.proj_0_data
        self.nproj0 = np.array(self.normproj0)
        self.normproj02 = self.nproj0/self.nproj0.max()*100
        
        return self.proj_0_scale, self.normproj02  
    
    
    def proj_0_hz(self, start, stop):
		
        self.proj_0_scale = self.uc0.hz_scale()
        self.proj_0_data = np.arange(0,len(self.data[0:,0]),1)
       
        for k in self.proj_0_data:
            test = sum(x for x in self.data[k,self.uc1(str(start)+' ppm'):self.uc1(str(stop)+' ppm')] if x > 0)
            self.proj_0_data[k]=test/10000
		
        return self.proj_0_scale, self.proj_0_data
    
    
    def proj_0_hzn(self, start, stop):
		
        self.proj_0_scale = self.uc0.hz_scale()
        self.proj_0_data = np.arange(0,len(self.data[0:,0]),1)
		
        self.ndata = self.data/self.data.max()*100
		
        for k in self.proj_0_data:
            test = sum(x for x in self.ndata[k,self.uc1(str(start)+' ppm'):self.uc1(str(stop)+' ppm')] if x > 0)
            self.proj_0_data[k]=test
            
        return self.proj_0_scale, self.proj_0_data
		
    
    def build2d(self,size=(5,4),dpi=100,color='k',linewidth=0.5,limits=((10,-2),(20,-4)),use_cmap=False,cmap='jet',vmin=0,vmax=1e7):
        
        self.xaxis, self.yaxis, self.int = self.plot2d()
        
        self.fig = plt.figure(dpi=dpi, figsize=size)
                
        self.proj_x_scale, self.proj_x_data = self.proj_1(start=limits[1][0],stop=limits[1][1])
        self.proj_y_scale, self.proj_y_data = self.proj_0(start=limits[0][0],stop=limits[0][1])
        
        self.ax1 = plt.subplot2grid((8, 10), (1, 1), colspan=9, rowspan=7)
        if use_cmap == False:
            plt.contour(self.xaxis, self.yaxis, self.int, self.c_levels,colors=color,linewidths=linewidth)
        else:
            plt.contour(self.xaxis, self.yaxis, self.int, self.c_levels,cmap=cmap,vmin=vmin,vmax=vmax,linewidths=linewidth)
        plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.minorticks_on()
        plt.xlim(limits[0])
        plt.ylim(limits[1])        
        
        self.ax2= plt.subplot2grid((8, 10), (0, 1), colspan=9, rowspan=1)
        plt.plot(self.proj_x_scale,self.proj_x_data,'-',color='k',linewidth=1)
        
        plt.xlim(limits[0][:])
        plt.axis('off')
                
        self.ax3 = plt.subplot2grid((8, 10), (1, 0), colspan=1, rowspan=7)
        plt.plot(-1.0*self.proj_y_data,self.proj_y_scale,'-',color='k',linewidth=1)
        
        plt.ylim(limits[1][:])
        plt.axis('off')
        
        return self.fig, self.ax1, self.ax2, self.ax3

        
    def build2d2(self,size=(5,4),dpi=100,color='k',linewidth=0.5,limits=((10,-2),(20,-4)),use_cmap=False,cmap='jet',vmin=0,vmax=100):
        
        self.xaxis, self.yaxis, self.int = self.plot2d()
        
        self.fig = plt.figure(dpi=dpi, figsize=size)
        
        self.normdata = self.data
        self.ndata = np.array(self.normdata)
        self.normdata2 = self.ndata/self.ndata.max()*100
        
        
        self.proj_x_scale, self.proj_x_data = self.proj_1_2(start=limits[1][0],stop=limits[1][1])
        self.proj_y_scale, self.proj_y_data = self.proj_0_2(start=limits[0][0],stop=limits[0][1])
        
        self.ax1 = plt.subplot2grid((8, 10), (1, 1), colspan=9, rowspan=7)
        if use_cmap == False:
            plt.contour(self.xaxis, self.yaxis, self.normdata2, self.c_levels,colors=color,linewidths=linewidth)
        else:
            plt.contour(self.xaxis, self.yaxis, self.normdata2, self.c_levels,cmap=cmap,vmin=vmin,vmax=vmax,linewidths=linewidth)
        plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.minorticks_on()
        plt.xlim(limits[0])
        plt.ylim(limits[1])        
        
        self.ax2= plt.subplot2grid((8, 10), (0, 1), colspan=9, rowspan=1)
        plt.plot(self.proj_x_scale,self.proj_x_data,'-',color='k',linewidth=1)
        
        plt.xlim(limits[0][:])
        plt.axis('off')
                
        self.ax3 = plt.subplot2grid((8, 10), (1, 0), colspan=1, rowspan=7)
        plt.plot(-1.0*self.proj_y_data,self.proj_y_scale,'-',color='k',linewidth=1)
        
        plt.ylim(limits[1][:])
        plt.axis('off')
        
        return self.fig, self.ax1, self.ax2, self.ax3        

    
    def buildrocsa(self,size=(5,4),dpi=100,color='k',linewidth=0.5,limits=((10,-2),(20,-4)),use_cmap=False,cmap='jet',vmin=0,vmax=100):
        
        self.xaxis, self.yaxis, self.int = self.plot_rocsa()
        
        self.fig = plt.figure(dpi=dpi, figsize=size)
        
        self.normdata = self.data
        self.ndata = np.array(self.normdata)
        self.normdata2 = self.ndata/self.ndata.max()*100

        
        self.proj_x_scale, self.proj_x_data = self.proj_1_rocsa_norm(start=limits[1][0],stop=limits[1][1])
        self.proj_y_scale, self.proj_y_data = self.proj_0_hzn(start=limits[0][0],stop=limits[0][1])
        
        self.ax1 = plt.subplot2grid((8, 10), (1, 1), colspan=9, rowspan=7)
        if use_cmap == False:
            plt.contour(self.xaxis, self.yaxis, self.normdata2, self.c_levels,colors=color,linewidths=linewidth)
        else:
            plt.contour(self.xaxis, self.yaxis, self.normdata2, self.c_levels,cmap=cmap,vmin=vmin,vmax=vmax,linewidths=linewidth)
        plt.tick_params(axis='x', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.tick_params(axis='y', which='both', labelleft=False, left=False, labelright=True, right=True,labelsize=10)
        plt.minorticks_on()
        plt.xlim(limits[0])
        plt.ylim(limits[1])        
        
        self.ax2= plt.subplot2grid((8, 10), (0, 1), colspan=9, rowspan=1)
        plt.plot(self.proj_x_scale,self.proj_x_data,'-',color='k',linewidth=1)
        
        plt.xlim(limits[0][:])
        plt.axis('off')
                
        self.ax3 = plt.subplot2grid((8, 10), (1, 0), colspan=1, rowspan=7)
        plt.plot(-1.0*self.proj_y_data,self.proj_y_scale,'-',color='k',linewidth=1)
        
        plt.ylim(limits[1][:])
        plt.axis('off')
        
        return self.fig, self.ax1, self.ax2, self.ax3

        
