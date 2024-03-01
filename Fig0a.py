import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf
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






rad=np.pi/180.

if True:

    Vels=-40*np.cos(np.arange(0,90,0.1)*rad)

    Acero=abs(Vels)<1
    Amas=abs(Vels)>1
    wvl=np.arange(-50,25,.1) 
    Dopp=6.
    perf=np.zeros(len(wvl))
    for i,w in enumerate(wvl):
        #v(z)= alfa z
        # x1=w
        # x2=w-Vels[Amas] 
        # tau1=-np.sqrt(np.pi)*Dopp/(2.*Vels[Amas])*(erf(x2/Dopp)-erf(x1/Dopp))
        # Tau=tau1    

        #v(z)=alfa (1-z)
        x1=w-Vels[Amas]  
        x2=w
        tau1=np.sqrt(np.pi)*Dopp/(2.*Vels[Amas])*(erf(x2/Dopp)-erf(x1/Dopp))
        Tau=tau1

        #v(z)=alfa* sqrt(z)
        # x1=w 
        # x2=w-Vels[Amas]   
        # tau1=-(2./Vels[Amas]**2)*(np.sqrt(np.pi)/2.)*w*Dopp*(erf(x2/Dopp)-erf(x1/Dopp))    
        # tau2=-(2./Vels[Amas]**2)*(Dopp**2/2.)*(np.exp(-x2**2/Dopp**2)-np.exp(-x1**2/Dopp**2))                
        # Tau=(tau1+tau2)

        #v(z)=alfa * sqrt(1-beta z)
        # beta=.2
        # x1=w-Vels[Amas]
        # x2=w-Vels[Amas]*np.sqrt(1.-beta)
        # tau1=(np.sqrt(np.pi)/(beta*Vels[Amas]**2))*w*Dopp*(erf(x2/Dopp)-erf(x1/Dopp))  
        # tau2=(Dopp**2/(beta*Vels[Amas]**2))*(np.exp(-x2**2/Dopp**2)-np.exp(-x1**2/Dopp**2))  
        # Tau=tau1+tau2 


        perf[i]=np.sum(len(Vels[Amas])*(np.exp(-Tau)))+np.sum(len(Vels[Acero])*(1.-np.exp(-(w-Vels[Acero])**2/Dopp**2)))

    #plt.subplot(1,2,2)
    plt.plot(wvl,perf/np.amax(perf),color='b',linewidth=2)
    plt.axvline(0,linestyle='--')
    plt.axvline(-40,linestyle='--')

    Vels=-40
    x1=wvl-Vels 
    x2=wvl
    tau1=np.sqrt(np.pi)*Dopp/(2.*Vels)*(erf(x2/Dopp)-erf(x1/Dopp))
    Tau=tau1
    plt.plot(wvl,np.exp(-Tau),'r--')

    Vels=-20
    x1=wvl-Vels 
    x2=wvl
    tau1=np.sqrt(np.pi)*Dopp/(2.*Vels)*(erf(x2/Dopp)-erf(x1/Dopp))
    Tau=tau1
    plt.plot(wvl,np.exp(-Tau),'r--')

    Vels=-10
    x1=wvl-Vels 
    x2=wvl
    tau1=np.sqrt(np.pi)*Dopp/(2.*Vels)*(erf(x2/Dopp)-erf(x1/Dopp))
    Tau=tau1
    plt.plot(wvl,np.exp(-Tau),'r--',label=r'$\theta=75,5^{\circ}$')
    # bipos,bival=bisector(wvl,perf)
    # plt.plot(bipos,bival)

plt.xlabel('Velocity (km/s)') 
plt.ylabel('Intensity') 
plt.show()
