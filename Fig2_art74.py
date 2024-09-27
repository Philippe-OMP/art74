import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.special import erf,roots_legendre
from scipy.optimize import fsolve

def bisector(wvl,perf):

    a1=np.amin(perf)
    abajo=np.argmin(perf)
    a2=np.amax(perf)
    margen=0.05*(a2-a1)
    bival=np.linspace(a1+margen,a2-margen,20)
    bipos=np.zeros(20)
    
    for i,x in enumerate(bival):
        este=np.argmin(abs(perf[0:abajo]-x))
        if perf[este]-x>0:
            otro=este+1
        else:
            otro=este-1
        p1=np.interp(x,[perf[este],perf[otro]],[wvl[este],wvl[otro]])

        este=abajo+np.argmin(abs(perf[abajo:]-x))
        if perf[este]-x>0:
            otro=este+1
        else:
            otro=este-1
        p2=np.interp(x,[perf[este],perf[otro]],[wvl[este],wvl[otro]])
        bipos[i]=0.5*(p2+p1)
        
    
    return bipos,bival






# Parabolic deceleration
def RT(wvl,vcos,Dopp,beta):
    x= lambda z: np.sqrt(1.-beta*z)
 
    if abs(vcos)>1:               
        tau1=np.sqrt(np.pi)*wvl*Dopp*(erf((wvl-vcos*x(1))/Dopp)-erf((wvl-vcos*x(0))/Dopp))  
        tau2=Dopp**2*(np.exp(-(wvl-vcos*x(1))**2/Dopp**2)-np.exp(-(wvl-vcos*x(0))**2/Dopp**2))  
        Tau=(tau1+tau2)/(beta*vcos**2)  
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-Tau)
    return perfil

def RT_linear(wvl,vcos,Dopp,beta):
    if abs(vcos)>1:
        tau1=np.sqrt(np.pi)*Dopp*(erf((wvl-vcos*beta)/Dopp)-erf((wvl/Dopp)))    
        Tau=tau1/(beta*2)
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-Tau)
    return perfil

def RT_constant(wvl,vcos,Dopp,beta):
    
    Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-Tau)
    return perfil


rad=np.pi/180.


if True: # Fig2a
    figname="Fig2a.png"
    transfer=RT_constant  # or RT_linear
else: # Fig2b
    figname="Fig2b.png"
    transfer=RT

V0=-40

wvl=np.arange(-50,25,.1) 
Dopp=6.
perf=np.zeros(len(wvl))

beta=1

Vels=-40
perf=transfer(wvl,Vels,Dopp,beta)   #cos mu to project onto disk
plt.plot(wvl,perf/perf.max(),'r--')

Vels=-40*np.cos(60*rad)
perf=transfer(wvl,Vels,Dopp,beta)   #cos mu to project onto disk
plt.plot(wvl,perf/perf.max(),'r--')

Vels=-40*np.cos(85*rad)
perf=transfer(wvl,Vels,Dopp,beta)  #cos mu to project onto disk
plt.plot(wvl,perf/perf.max(),'r--')


for mu in np.arange(0,90,0.1)*rad:
    perf=perf+transfer(wvl,V0*np.cos(mu),Dopp,beta)*np.cos(mu)   #cos mu to project onto disk

#plt.subplot(1,2,2)
plt.plot(wvl,perf/np.amax(perf),color='b',linewidth=2)
plt.axvline(0,linestyle='--')
plt.axvline(-40,linestyle='--')




plt.xlabel('Velocity (km/s)') 
plt.ylabel('Intensity') 
plt.savefig(figname)
plt.show()

