#-*-coding:utf-8-*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
#from mpl_toolkits.basemap import Basemap 

class Ylm():

	def __init__(self,lmax,star=''):
		self.lmax=lmax


		self.longitude = np.arange(0,361,2)
		self.latitude = np.arange(-10, 90, 1)	#First number less than 0....peek on the other side
		#self.latitude = np.arange(0, 90, 1)	 #No peeking
		self.mu=(90 - self.latitude)*np.pi/180 
		self.teta=(self.longitude*np.pi/180) 

		self.tt,self.mm=np.meshgrid(self.teta,self.mu)
	# 
	# r=np.linspace(0,1,200)
	# theta=np.linspace(0,2*np.pi,1800)
	# rr,th=np.meshgrid(r,theta)
	# x=np.fix(200*rr*np.cos(th)+220)
	# y=np.fix(200*rr*np.sin(th)+220)

		file="/Users/art2/Betelgeuse/avQU/n20131127_q.av"
		self.V0=40.	 #40 for Betelgeuse; 50 for muCeph? ; 30 for CETau
		self.Vstar=27  #km/s Proper velocity of Betelgeuse(25) ; 23 muCeph, 23.75 CE Tau
		if star=="CETau":
			dir='/Users/arturo/CETau/profs/'
			file=dir+'n20150310_q.s' #CE Tau
			self.V0=30.	 #40 for Betelgeuse; 50 for muCeph? ; 30 for CETau
			self.Vstar=23.75  #km/s Proper velocity of Betelgeuse(25) ; 23 muCeph, 23.75 CE Tau
		if star=="muCeph":
			dir='/Users/arturo/muCeph/profiles/' #'../avQU/'
			file=dir+'n20150711_q.s' #mu Ceph
			self.V0=50.	 #40 for Betelgeuse; 50 for muCeph? ; 30 for CETau
			self.Vstar=23  #km/s Proper velocity of Betelgeuse(25) ; 23 muCeph, 23.75 CE Tau
		if star=="SUPer" or star=="SU_Persei":
			dir='/Users/arturo/SUPersei/profiles/' #'../avQU/'
			file=dir+'n20161210_q.s' #SUPer
			self.V0=30.	 #40 for Betelgeuse; 50 for muCeph? ; 30 for CETau
			self.Vstar=-40	#km/s Proper velocity of Betelgeuse(25) ; 23 muCeph, 23.75 CE Tau
		if star=="AZCyg":
			dir='/Users/arturo/AZCygni/profiles/' #'../avQU/'
			file=dir+'n20180813_q.s' #SUPer
			self.V0=30.	 #40 for Betelgeuse; 50 for muCeph? ; 30 for CETau
			self.Vstar=-5  #km/s Proper velocity of Betelgeuse(25) ; 23 muCeph, 23.75 CE Tau
		
	#files=sorted([f for f in os.listdir('.') if f.endswith('q.av')])
		Stokes,dimx,dimy=self.read_Stokes(file)
		
		self.wvl=Stokes[:,0]
		self.fwhm=10.  #Avg standard deviation of observed profiles in Betelgeuse. But.....
		
		self.cuantos=lmax*(lmax+2)+1
	# vbase=v0*np.cos(rr*np.pi/2.)
	# PA=(np.pi/2.+tt).flatten()
	# Raybase=rr*(1.-np.cos(rr*np.pi/2.)**2)

	def Yreal(self,PL,m):
		if m>0:
			Yr=np.sqrt(2)*(-1)**m*PL.real
		elif m<0:
			Yr=np.sqrt(2)*(-1)**m*PL.imag
		else:
			Yr=PL.real
		return Yr
		
	def One_Ylm(self,cual):
	
		l=self.noll_n(cual)
		m=self.noll_m(cual)
		PL= sph_harm(abs(m),l,self.tt,self.mm)
		One=self.Yreal(PL,m)
		return One

	def Create_Ylm_star(self,coeff):

		L=np.arange(1,self.lmax+1,1) #Add +4  to add high-freq noise

		SUM=0.
		
		sigmafactor=[0,0,0,0.01,0.1,0.1,0.1,0,2,10,10]

		for l in L:

			for m in range(0,l+1,1):

				#coeff=uniform(-1,1,size=2) 
				#coeff=coeff/np.sum(abs(coeff)) 
				#print(l,m)
				PL= sph_harm(abs(m),l,self.tt,self.mm)
				cual=self.noll(l,m)
				
				
				
				if l>self.lmax:
					cc=0.#np.random.randn()*np.mean(coeff[10:])
				else:
					cc=coeff[cual]
				#cc=cc*(1.+sigmafactor[l]*np.random.randn())
				SUM=SUM+cc*self.Yreal(PL,m)#coeff[1,cual]*PL.imag
				
				if m>0:
					cual=self.noll(l,-m)
					cc=coeff[cual] ####???????    Esto no estaba cuando he corregido en XII-19. Por que?
					SUM=SUM+cc*self.Yreal(PL,-m)#coeff[1,cual]*PL.imag
				
				#SUM.append(hs) 
				#print 'm=',m,'-----','coeff=',coeff,'\n'
		 


		B=SUM#np.abs(SUM)
		return B
		
	def Ylm_Vels(self,Star):

		m0=np.amin(abs(Star))
		dif=np.ptp(abs(Star))
		V=-self.V0*np.cos(self.mm)*(-1+2.*np.tanh(2.1972*(abs(Star)-m0)/dif)) # The darkest 25% goes down, the rest up
		return V
		
	def Compute_QU(self,Star):


		c=self.fwhm/(2*np.sqrt(np.log(2)))


		perf=np.zeros(len(self.wvl))
		perfQ=np.zeros(len(self.wvl))
		perfU=np.zeros(len(self.wvl))
	#	for (i,j),theta in np.ndenumerate(tt):
	#		V=V_0*np.cos(mm[i,j])
	#		perf=perf+(1.-B[i,j]*np.exp(-(wvl-V)**2/(c**2)))*np.sin(mm[i,j])
	#		perfQ=perfQ+B[i,j]*np.sin(mm[i,j])**2*np.cos(2*theta)*np.exp(-(wvl-V)**2/c**2)*np.sin(mm[i,j])
	#		perfU=perfU+B[i,j]*np.sin(mm[i,j])**2*np.sin(2*theta)*np.exp(-(wvl-V)**2/c**2)*np.sin(mm[i,j])
	#	
	
		#vbase=-self.V0*np.cos(self.mm)*abs(Star)/np.amax(abs(Star))   # brightness only. Peek on the other side
		m0=np.amin(abs(Star))
		dif=np.ptp(abs(Star))
		#vbase=-self.V0*np.cos(self.mm)*(-1+2.*(abs(Star)-m0)/dif)	# Dark goes down, Bright goes up. Linear dependence 
		#vbase=-self.V0*np.cos(self.mm)*(-2+3.*np.tanh(5.493*(abs(Star)-m0)/dif)) # The darkest 10% goes down, the rest up
		#vbase=-self.V0*np.cos(self.mm)*(-1+2.*np.tanh(2.7465*(abs(Star)-m0)/dif)) # The darkest 20% goes down, the rest up
		vbase=-self.V0*np.cos(self.mm)*(-1+2.*np.tanh(2.1972*(abs(Star)-m0)/dif)) # The darkest 25% goes down, the rest up
		#vbase=-self.V0*np.cos(self.mm)*(-1+2.*np.tanh(1.6479*(abs(Star)-m0)/dif)) # The darkest 33% goes down, the rest up
		#Towards us is negative and bright....
		
		

		RayQ=abs(Star)*np.sin(self.mm)**2*np.cos(2*self.tt+np.pi)#*np.sin(self.mm) #sin**2 mu Rayleigh rdr=sin	mu integration
		RayU=abs(Star)*np.sin(self.mm)**2*np.sin(2*self.tt+np.pi)#*np.pi is to correct respect to Narval polarimeter (see Auriere 2016)
		for i,w in enumerate(self.wvl):
			PSF=np.exp(-(w-vbase-self.Vstar)**2/(c**2))
			perf[i]=np.sum(abs(Star)*(1.-0.8*PSF*np.sin(self.mm))/np.amax(abs(Star)))	# rdr=sin integration
			perfQ[i]=np.sum(RayQ*PSF)
			perfU[i]=np.sum(RayU*PSF)

		
		perfQ=perfQ/perf.max()
		perfU=perfU/perf.max()
		perf=perf/perf.max()
		return perfQ,perfU,perf


	def plot_Ylm(self,coef):
	
		Star=self.Create_Ylm_star(coef)
		pQ,pU,pI=self.Compute_QU(Star)
		
		plt.subplot(1,2,1)
		map= Basemap(projection='ortho',lat_0=90, lon_0=0)
		longs,lats=np.meshgrid(self.longitude,self.latitude)
		#longg,latt=map(longs,lats)
		clevs=np.linspace(0,np.amax(Star),100)
		map.contourf(longs,lats,Star,clevs,latlon=True)
		#map.colorbar()
		
		ax=plt.subplot(1,2,2)
		ax.plot(self.wvl,pQ,lw=2)
		ax.plot(self.wvl,pU,lw=2)
		ax1=ax.twinx()
		ax1.plot(self.wvl,pI,color='black')
		plt.xlim((-100,100))
		plt.draw()
		return

	def comparing_plot(self,coef1,coef2):
	
		
		
		Star=self.Create_Ylm_star(coef1)
		pQ1,pU1,pI=self.Compute_QU(Star)
		
		plt.subplot(2,2,3)
		map= Basemap(projection='ortho',lat_0=90, lon_0=0)
		longs,lats=np.meshgrid(self.longitude,self.latitude)
		#longg,latt=map(longs,lats)
		a1=np.amax(Star)
		a0=np.amin(Star)
		clevs=np.linspace(a0,a1,100)
		map.contourf(longs,lats,Star,clevs,latlon=True)
		
		Star=self.Create_Ylm_star(coef2)
		pQ2,pU2,pI=self.Compute_QU(Star)
		
		plt.subplot(2,2,4)
		map= Basemap(projection='ortho',lat_0=90, lon_0=0)
		clevs=np.linspace(a0,max([a1,np.amax(Star)]),100)
		map.contourf(longs,lats,Star,clevs,latlon=True)
		
		plt.subplot(2,2,1)

		
		plt.plot(self.wvl,pQ1)
		plt.plot(self.wvl,pQ2,linewidth=2)
		plt.xlim((-100,100))
		
		ax=plt.subplot(2,2,2)
		plt.plot(self.wvl,pU1)
		plt.plot(self.wvl,pU2,linewidth=2)
		plt.xlim((-100,100))
		
		
		
		plt.draw()
		return

	def read_Stokes(self,fich):
		f=open(fich,'r')
		f.readline()
		linea=f.readline()
		linea=linea.strip()
		dimx,dimy=map(int,linea.split())

		resto=f.read()
		DatosQ=np.asarray(list(map(float,resto.split())))
		DatosQ=DatosQ.reshape((dimx,dimy+1))
		f.close()

		ind=fich.index('_q')
		s=list(fich)
		s[ind+1]='u'
		file=''.join(s)

		f=open(file,'r')
		f.readline()
		linea=f.readline()
		linea=linea.strip()
		dimx,dimy=map(int,linea.split())

		resto=f.read()
		DatosU=np.asarray(list(map(float,resto.split())))
		DatosU=DatosU.reshape((dimx,dimy+1))
		f.close()
		Stokes=np.zeros((dimx,6))
		Stokes[:,0]=DatosQ[:,0]	 #Wavelength
		Stokes[:,1]=DatosQ[:,1]	 #Intensity from Q
		Stokes[:,2]=DatosQ[:,3]	 #Q
		Stokes[:,3]=-DatosU[:,3]  # U
		Stokes[:,4]=DatosQ[:,5] # Null in Q
		Stokes[:,5]=DatosU[:,5] # Null in U
		return Stokes,dimx,dimy
	

	
	def noll(self,n,m):
			idx=0
		
			if np.abs(m)<=n:
				idx=(n**2-1) +1 #This is for Yl m=0
				idx=idx+2*np.abs(m)
				if m>0:
					idx=idx-1
			
			return int(idx)

	def noll_n(self,j):

			idx=np.floor(np.sqrt(j))
			
			return int(idx)

	def noll_m(self,j):


			n=self.noll_n(j)
			idxabs=np.ceil((j-self.noll(n,0))/2.)
			m=idxabs
			if (m==(j-self.noll(n,0))/2.):
				m=-m
			return int(m)
	def Ambiguity(self,Y2):
		for j in range(self.cuantos):
			m=self.noll_m(j)
			Y2[j]=Y2[j]*(-1)**m
		return Y2

