### Version 06.07.2017
### works on ### Python 3.6.1 -- Matplotlib 2.0.0 -- Numpy 1.12.1 -- Scipy 0.19.0
################################
### Alexander Riss & Mathias Pörtner : Data loading and header parsing
### Mathias Pörtner : Adding graph abilities and plotting spectral line 'maps', GUI
### Domenik Zimmermann	: Documentation and dependency testing
################################
### For dependencies on Linux Mint 18 Cinnamon 64-bit
################################
### For dependencies on Windows
### Install Anaconda (tested for Anaconda3)
### Run an Ipython console with this file
################################
# Copy the script to directory of the data files
# Edit the topographic image you want to show with gwyddion and save it as ASCII data matrix (.txt) with the option "Add informational header" 
# Important: Do not trim images, and record them in the same size as u did the spectra on, otherwise thepositions will not match
################################
### Chose image color and spectra color ###
contrast_spec='afmhot'
contrast_topo='viridis'
### Choose scalebar length in nm ###
scalebar_length=5
### Average # points of spectra ###
average_specs=5
################################

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pylab as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import re
import glob
import numpy as np
import tkinter as tk
from tkinter import ttk
from scipy.optimize import curve_fit as cv

fontna=16
fontnu=12

totalfig = Figure(figsize=(12,6), dpi=100, tight_layout=True)
aplot = totalfig.add_subplot(121)
bplot = totalfig.add_subplot(122)

class SeaofBTCapp(tk.Tk):
	def __init__(self, *args, **kwargs):
		tk.Tk.__init__(self, *args, **kwargs)
		
		container = tk.Frame(self)
		container.pack(side='top', fill='both', expand=True)
		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)
		
		
		self.frames={}
	
		frame = PageThree(container,self)
		self.wm_title("Plot line Spectra")
		self.frames[PageThree] = frame
		frame.grid(row=0,column=0,sticky='nsew')
		
     
class PageThree(tk.Frame):
	def __init__(self,parent,controller):
		#Initialize GUI for choosing files
		tk.Frame.__init__(self,parent)
					
		self.widget=None
		self.toolbar=None
		
		#Choose topo file (.txt)
		labe1=tk.Label(self, text='Image-file').grid(row=0,column=0,sticky='w')
		choicesIma = glob.glob('*.txt')
		choicesIma.sort()
		self.variableIma = tk.StringVar(self)
		self.variableIma.set(choicesIma[0])
		wIma = ttk.Combobox(self, textvariable=self.variableIma, values=choicesIma, width=25)
		wIma.grid(row=0,column=1,sticky='w')
		
		#Choose file of spectra along line (.L*.VERT) (only first spectrum is displayed
		labe2=tk.Label(self, text='Spectra-file').grid(row=1,column=0,sticky='w')
		choicesSpec = glob.glob('*L0001.VERT')
		choicesSpec.sort()
		self.variableSpec = tk.StringVar(self)
		self.variableSpec.set(choicesSpec[0])
		wSpec = ttk.Combobox(self, textvariable=self.variableSpec, values=choicesSpec, width=25)
		wSpec.grid(row=1,column=1,sticky='w')
		
		#Load data button
		loadbuttonIma=tk.Button(self,text='Load',command=self.plotimage, width=12)
		loadbuttonIma.grid(row=0,column=2)
		
		#Quit button
		button = tk.Button(self, text='Quit', command=self._quit, width=12)
		button.grid(row=0,column=3)
		
		self.imau = tk.DoubleVar(self)
		self.imao = tk.DoubleVar(self)
		self.spu = tk.DoubleVar(self)
		self.spo = tk.DoubleVar(self)
	
	def _quit(self):
		self.quit()
		self.destroy()

	def string_simplify(self,str):
		#simplifies a string (i.e. removes replaces space for "_", and makes it lowercase
		return str.replace(' ','_').lower()

	def laden_spec(self,data):
		#Reads .VERT file for spectral information, returns data-array for Voltage, dIdV and spectrum position	
		headersp = {}
		f = open(data, encoding='utf-8', errors='ignore')
		headersp_ended = False
		caption = re.compile(':*:')
		key = ''
		contents = ''
		while not headersp_ended:
			line = f.readline()
			if not line: break
			if line[0:4] == "    ":
				parts=line.split()
				posi=np.array([float(parts[-2]),float(parts[-1])],float)
				headersp_ended = True
			else:
				parts = line.split('=')
				if len(parts)!=2: continue
				key, contents = parts
				line = line.strip()
				key = self.string_simplify(key)
				headersp[key] = contents.strip()
		f.close()
		
		global pixelsize
		dacstep=np.array([float(headersp['delta_x_/_delta_x_[dac]']),float(headersp['delta_y_/_delta_y_[dac]'])],float)
		pixelsize=np.array([float(headersp['num.x_/_num.x']),float(headersp['num.y_/_num.y'])],float)
		imagesize=np.array([float(headersp['length_x[a]']),float(headersp['length_y[a]'])],float)
		
		posi=posi/dacstep
		posi[0]=(pixelsize[0]/2.0+posi[0])*imagesize[0]/pixelsize[0]/10
		posi[1]=(pixelsize[1]-posi[1])*imagesize[1]/pixelsize[1]/10
		
		
		
		A=np.genfromtxt(data,delimiter='	',skip_header=212,skip_footer=0)
		U=A[:,3]
		dIdU=A[:,2]
		return(U,dIdU,posi)

	def laden_image(self,data):
		#Reads .txt file from Gwyddion for topographic information, returns 2D data-array for topography and the size of the image in nm
		header = {}
		f = open(data, encoding='utf-8', errors='ignore')
		header_ended = False
		caption = re.compile(':*:')
		key = ''
		contents = ''
		while not header_ended:
			line = f.readline()
			if not line: break
			if line[0] != "#":
				header_ended = True
			else:
				parts = line.split(':')
				if len(parts)!=2: continue
				key, contents = parts
				line = line.strip()
				key = self.string_simplify(key[2:])
				header[key] = contents[:-4].strip()
		f.close()
		
		ext=np.array([float(header['width']),float(header['height'])],float)
		
		X=np.loadtxt(data)*1e10
		
		return(X,ext)

	def minmax(self,X):
		#calculation of the minimal and maximal values of the 2D arrays (for setting the contrast)
		mini=100.0
		maxi=0.0
		for i in X:
			if max(i)>maxi:
				maxi=max(i)
			if min(i)<mini:
				mini=min(i)
		mean=np.mean([mini,maxi])
		return(mini,mean,maxi)

	def update_imau(self,imau):
		global canvas
		global imagepl
		vuntenima = self.imau.get()
		imagepl.set_clim(vmin=vuntenima)
		canvas.show()
	
	def update_imao(self,imao):
		global canvas
		global imagepl
		vobenima = self.imao.get()
		imagepl.set_clim(vmax=vobenima)
		canvas.show()
	
	def update_spu(self,spu):
		global canvas
		global specpl
		vuntensp = self.spu.get()
		specpl.set_clim(vmin=vuntensp)
		canvas.show()
		
	def update_spo(self,spo):
		global canvas
		global specpl
		vobensp = self.spo.get()
		specpl.set_clim(vmax=vobensp)
		canvas.show()
		
	def reset(self):
		imagepl.set_clim(vmin=self.imau0,vmax=self.imao0)
		specpl.set_clim(vmin=self.spu0,vmax=self.spo0)
		canvas.show()
		
	def sel(self):
		untenima=self.imau.get()
		obenima=self.imao.get()
		untensp=self.spu.get()
		obensp=self.spo.get()
		np.savetxt(self.data_name[:-4]+'.csv',[untenima,obenima,untensp,obensp],delimiter=',')
		
	def saveima(self):
		totalfig.savefig(self.data_name[:-4]+'.pdf')
		
	def ocon(self):
		con=np.loadtxt(self.data_name[:-4]+'.csv',delimiter=',')
		imagepl.set_clim(vmin=con[0],vmax=con[1])
		specpl.set_clim(vmin=con[2],vmax=con[3])
		canvas.show()
		
	def averagespec(self,matrixx,matrixy,ave):
		#average spectra
		matrixyneu=[]
		matrixxneu=[]
		for n,i in enumerate(matrixy):
			matrixyneu.append([])
			matrixxneu.append([])
			for m,j in enumerate(i):
				if m%ave==(ave-1):
					matrixyneu[-1].append(sum(matrixy[n][m-(ave-1):m])/ave)
					matrixxneu[-1].append(sum(matrixx[n][m-(ave-1):m])/ave)
		matrixyneu=np.array(matrixyneu,float)
		matrixxneu=np.array(matrixxneu,float)
		return(matrixxneu,matrixyneu)
		
	def plotimage(self):
		#Plot image and set up final GUI
		global canvas
		global unteima
		global obenima
		global untensp
		global obensp
		global imagepl
		global specpl
		
		s=self.variableSpec.get()
		spec=glob.glob(s[:-9]+'*.VERT')
		spec.sort()
		self.data_name=self.variableIma.get()
		
		ima,imagesize=self.laden_image(self.data_name)
		
		ext=[]
		matrixx=[]
		matrixy=[]
		spec_posi=[]
		for i in spec:
			x,y,posi=self.laden_spec(i)
			matrixy.append(y)
			matrixx.append(x)
			spec_posi.append(posi)
			if i==spec[0] or i==spec[-1]:
				ext.append(posi)
		matrixx=np.array(matrixx,float)
		matrixy=np.array(matrixy,float)

		line_length=np.sqrt((ext[0][0]-ext[1][0])**2+(ext[0][1]-ext[1][1])**2)
		
		#normalize specs
		maxi=[]
		globmaxy=0
		for n,i in enumerate(matrixy):
			maxi.append(0)
			for m,j in enumerate(i):
				if matrixx[n][m]>-700 and matrixx[n][m]<0 and matrixy[n][m]>maxi[n]:
					maxi[n]=matrixy[n][m]
			if max(i)>globmaxy:
				globmaxy=max(i)
				
		matrixxn,matrixyn=self.averagespec(matrixx,matrixy,average_specs)
	
		if self.widget:
			self.widget.destroy()
		if self.toolbar:
			self.toolbar.destroy()
		aplot.clear()
		bplot.clear()
		
		untenima,mitteima,obenima=self.minmax(ima)
		untensp,mittesp,obensp=self.minmax(matrixyn)
		
		self.imau0=untenima
		self.imao0=obenima
		self.spu0=untensp
		self.spo0=obensp
		
		#Plot Topography
		imagepl=aplot.imshow(ima,cmap=contrast_topo,extent=[0,imagesize[0],0,imagesize[1]],vmin=untenima,vmax=obenima)
		aplot.set_xticks([])
		aplot.set_yticks([])
		for pos in spec_posi:
			aplot.plot(pos[0],pos[1],'.',color='white',ms=1)
		for x in np.linspace(1,1+scalebar_length,1000):
			aplot.plot(x,1,'.',color='white',ms=1)
		
		#Plot spectra map
		specpl=bplot.imshow(matrixy.T,cmap=contrast_spec,extent=[0,line_length,min(matrixx[0]),max(matrixx[0])],aspect='auto',vmin=untensp,vmax=obensp)
		bplot.set_xlabel('Distance x [nm]',fontsize=fontna)
		bplot.set_ylabel('Bias voltage [mV]',fontsize=fontna)
		bplot.set_yticks([-1000,0,1000])
		
		#Put Plot in GUI
		canvas=FigureCanvasTkAgg(totalfig, self)
		canvas.show()
		self.widget=canvas.get_tk_widget()
		self.widget.grid(row=2,columnspan=4)
		
		#Sliders for contrast
		sluima=tk.Scale(self, from_=untenima-(mitteima-untenima), to=obenima,resolution=0.01,variable=self.imau,command=self.update_imau,orient=tk.HORIZONTAL)
		sluima.set(untenima)
		sluima.grid(row=3,column=1)
		
		labeliu=tk.Label(self, text='Lower value image')
		labeliu.grid(row=3,column=0)
		
		sloima=tk.Scale(self, from_=untenima, to=obenima+(obenima-mitteima),resolution=0.01,variable=self.imao,command=self.update_imao,orient=tk.HORIZONTAL)
		sloima.set(obenima)
		sloima.grid(row=4,column=1)
		
		labelio=tk.Label(self, text='Upper value image')
		labelio.grid(row=4,column=0)
		
		slusp=tk.Scale(self, from_=untensp-(mittesp-untensp), to=obensp,resolution=0.01,variable=self.spu,command=self.update_spu,orient=tk.HORIZONTAL)
		slusp.set(untensp)
		slusp.grid(row=3,column=3)
		
		labelsu=tk.Label(self, text='Lower value spectra')
		labelsu.grid(row=3,column=2)
		
		slosp=tk.Scale(self, from_=untensp, to=obensp+(obensp-mittesp),resolution=0.,variable=self.spo,command=self.update_spo,orient=tk.HORIZONTAL)
		slosp.set(obensp)
		slosp.grid(row=4,column=3)
		
		labelso=tk.Label(self, text='Upper value spectra')
		labelso.grid(row=4,column=2)
		
		#Save contrast button
		buttonsav = tk.Button(self, text="Save contrast", command=self.sel, width=12)
		buttonsav.grid(row=3,column=8)
		
		#Old contrast button
		buttonoc = tk.Button(self, text="Old contrast", command=self.ocon, width=12)
		buttonoc.grid(row=4,column=8)
		
		#Save image button
		buttonsavima = tk.Button(self, text="Save image", command=self.saveima, width=12)
		buttonsavima.grid(row=1,column=2)
		
		#Reset contrast button
		buttonres = tk.Button(self, text='Reset', command=self.reset, width=12)
		buttonres.grid(row=1,column=3)
		
app = SeaofBTCapp()
app.mainloop()
