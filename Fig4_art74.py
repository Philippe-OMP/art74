import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf
from scipy.optimize import fsolve



# Parabolic deceleration
def RT(wvl,vcos,Dopp,beta):
    x= lambda z: np.sqrt(1.-beta*z)
 
    if abs(vcos)>1:               
        tau1=np.sqrt(np.pi)*wvl*Dopp*(erf((wvl-vcos*x(1))/Dopp)-erf((wvl-vcos*x(0))/Dopp))  
        tau2=Dopp**2*(np.exp(-(wvl-vcos*x(1))**2/Dopp**2)-np.exp(-(wvl-vcos*x(0))**2/Dopp**2))  
        Tau=(tau1+tau2)/(beta*vcos**2)   # This vcos**2 is what justifies the separation of cases for vcos ~ 0
        
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-Tau)
    return perfil

rad=np.pi/180.
gauge=60  # Amplitude factor between observations and database	

Y=Ylm(5,star='Betelgeuse')
ficheroSolucion='Inversion_Y5.pkl'

fout=open(ficheroSolucion,'rb')
Solucion=pickle.load(fout,encoding="latin1")
fout.close()
#Solucion.append({'Date':fecha,'Ylm':Coef,'StokesQ':ProfilQ,'StokesU':ProfilU,'Data':Stokes,'Chi2':Chi2,'PCAdist':dist,'Gauge':gauge})	
cuantos=len(Solucion)

cual=2

S=Solucion[cual]
print(S['Date'],Solucion[10]['Date'])
Y0=S['Ylm']
Star=Ylm.Create_Ylm_star(Y,Y0)
Star=Star/np.sum(abs(Star))
Starmem=Ylm.Create_Ylm_star(Y,Y0)
Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-30,40)   #  LOS projected
longs,lats=np.meshgrid(Y.longitude,Y.latitude)  # in degrees

radio=np.where(longs==0)
Ecuador=(longs==0)

wvl=np.arange(-50,25,.1) 
Dopp=6.
perf=np.zeros(len(wvl))
beta=1
perfsuma=0.

for punto,long in enumerate(Y.longitude):
    if long==0:
        for plat,lat in enumerate(Y.latitude):
            if lat>0:

                perf=abs(Star[plat,punto])*RT(wvl,Vels[plat,punto],Dopp,beta)#*np.sin(lat*rad) # cos mu, latitude is 0 at limb
                plt.plot(wvl,perf/np.amax(perf),color='r',alpha=0.1)
                perfsuma=perfsuma+perf



plt.plot(wvl,perfsuma/perfsuma.max(),'b')
plt.axvline(0,linestyle='--')
plt.axvline(-40,linestyle='--')
plt.xlabel('Wavelength (km/s)')
plt.ylabel('Intensity')
plt.savefig('Fig4.png')

plt.show()
