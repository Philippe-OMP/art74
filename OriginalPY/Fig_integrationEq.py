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
    margen=0.01*(a2-a1)
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
gauge=60  # Amplitude factor between observations and database	

Y=Ylm(5,star='Betelgeuse')
ficheroSolucion='Inversion_Y5.pkl'

fout=open(ficheroSolucion,'rb')
Solucion=pickle.load(fout,encoding="latin1")
fout.close()
#Solucion.append({'Date':fecha,'Ylm':Coef,'StokesQ':ProfilQ,'StokesU':ProfilU,'Data':Stokes,'Chi2':Chi2,'PCAdist':dist,'Gauge':gauge})	
cuantos=len(Solucion)

cual=0

S=Solucion[cual]
print(S['Date'],Solucion[10]['Date'])
Y0=S['Ylm']
Star=Ylm.Create_Ylm_star(Y,Y0)
Star=Star/np.sum(abs(Star))
Starmem=Ylm.Create_Ylm_star(Y,Y0)
Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-30,30)
longs,lats=np.meshgrid(Y.longitude,Y.latitude)

radio=np.where(longs==0)
Ecuador=(longs==0)


# Star=Star[Ecuador]
# Vels=Vels[Ecuador]

Acero=abs(Vels)<1
Amas=abs(Vels)>1
wvl=np.arange(-50,25,.1) 
Dopp=6.
perf=np.zeros(len(wvl))

beta=1
perfsuma=0.
for punto,long in enumerate(Y.longitude):
    if long==0:
        for plat,lat in enumerate(Y.latitude):
            if lat>0:
                x1=wvl-Vels[plat,punto]
                x2=wvl-Vels[plat,punto]*np.sqrt(1.-beta)
                tau1=(np.sqrt(np.pi)/(beta*Vels[plat,punto]**2))*wvl*Dopp*(erf(x2/Dopp)-erf(x1/Dopp))  
                tau2=(Dopp**2/(beta*Vels[plat,punto]**2))*(np.exp(-x2**2/Dopp**2)-np.exp(-x1**2/Dopp**2))  
                Tau=tau1+tau2 
                if abs(Vels[plat,punto])<1:
                    perf=abs(Star[plat,punto])*(1.-0.4*np.exp(-(wvl-Vels[plat,punto])**2/Dopp**2))
                else:
                    perf=abs(Star[plat,punto])*(np.exp(-0.4*Tau))
                plt.plot(wvl,perf/np.amax(perf),color='r',alpha=0.1)
                perfsuma=perfsuma+perf


perf=perf/np.amax(perf)
plt.plot(wvl,perf+cual*0.01,'b')
plt.axvline(0,linestyle='--')
plt.axvline(-40,linestyle='--')
plt.xlabel('Wavelength (km/s)')
plt.ylabel('Intensity')
plt.savefig('Fig_integrationEquateur.png')

plt.show()
