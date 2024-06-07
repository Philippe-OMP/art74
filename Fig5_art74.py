import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf,roots_legendre
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






# Parabolic deceleration
def RT(wvl,vcos,Dopp,beta):
    x= lambda z: np.sqrt(1.-beta*z)
 
    if abs(vcos)>1:               
        tau1=np.sqrt(np.pi)*wvl*Dopp*(erf((wvl-vcos*x(1))/Dopp)-erf((wvl-vcos*x(0))/Dopp))  
        tau2=Dopp**2*(np.exp(-(wvl-vcos*x(1))**2/Dopp**2)-np.exp(-(wvl-vcos*x(0))**2/Dopp**2))  
        Tau=(tau1+tau2)/(beta*vcos**2)   # This vcos**2 is what justifies the separation of cases for vcos ~ 0
        
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-0.45*Tau)
    return perfil

def coord2long(theta,phi):
    rad=np.pi/180.
    lat=coord2lat(theta,phi)
    #sinlong=np.sin(theta)*np.sin(phi)/np.cos(lat)
    #return np.arcsin(sinlong)
    coslong=np.cos(theta)/np.cos(lat)
    if coslong>1:
        print(theta/rad,phi/rad,coslong)
        coslong=1.
    if coslong<-1:
        print(theta/rad,phi/rad,coslong)
        coslong=-1.
    return np.arccos(coslong)

def coord2lat(theta,phi):
    #sinlat=np.sin(phi)*np.sin(theta)  # phi (0,pi)
    sinlat=np.cos(phi)*np.sin(theta)  # phi (-pi/2,pi/2)
    return np.arcsin(sinlat)
rad=np.pi/180.
gauge=60  # Amplitude factor between observations and database	

Y=Ylm(5,star='Betelgeuse')
ficheroSolucion='Inversion_Y5.pkl'

fout=open(ficheroSolucion,'rb')
Solucion=pickle.load(fout,encoding="latin1")
fout.close()
#Solucion.append({'Date':fecha,'Ylm':Coef,'StokesQ':ProfilQ,'StokesU':ProfilU,'Data':Stokes,'Chi2':Chi2,'PCAdist':dist,'Gauge':gauge})	
cuantos=len(Solucion)

for cual in range(10):

    S=Solucion[cual]
    Y0=S['Ylm']
    Star=Ylm.Create_Ylm_star(Y,Y0)
    Star=Star/np.sum(abs(Star))
    Starmem=Ylm.Create_Ylm_star(Y,Y0)
    Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-30,20) #  LOS projected
    longs,lats=np.meshgrid(Y.longitude,Y.latitude)  # in degrees
    



    wvl=np.arange(-50,25,.1) 
    
    perf=np.zeros(len(wvl))
    if True:
        beta=.1 # 0.2 for RW Cep, 1 for Betelgeuse
        Dopp=10.
    else:
        beta=1.
        Dopp=6.
    perfsuma=0.

   
    nodos=20
    roots,weights=roots_legendre(nodos) 
    phis=np.pi*np.arange(nodos+1)/nodos-np.pi/2
    thetas=np.arccos(roots) 
    #print(phis/rad,thetas/rad)
    
    for phi in phis:
        for i,theta in enumerate(thetas):
            lat=coord2lat(theta,phi)
            estelat=np.argmin(abs(Y.latitude*rad-lat))
            az=coord2long(theta,phi)
            estelong=np.argmin(abs(Y.longitude*rad-az))
            #print(abs(Star[estelat,estelong]),Vels[estelat,estelong],np.sin(lat))
            #print(phi/rad,theta/rad,lat/rad,az/rad)
            perf=abs(Star[estelat,estelong])*RT(wvl,Vels[estelat,estelong],Dopp,beta)*np.sin(theta)*np.cos(phi)
            perfsuma=perfsuma+(np.pi/nodos)*weights[i]*perf

    
    # stop  

    # for punto,long in enumerate(Y.longitude):
        
    #     for plat,lat in enumerate(Y.latitude):
    #         if lat>0:
    #             cosmu=np.sin(lat*rad) # latitude is 0 at limb, and 90 at diskcenter
    #             perf=abs(Star[plat,punto])*RT(wvl,Vels[plat,punto],Dopp,beta)#*cosmu # cos mu for projection onto disk
    #             #plt.plot(wvl,perf/np.amax(perf),color='r',alpha=0.1)
    #             perfsuma=perfsuma+perf

    plt.subplot(1,2,1)
    perfsuma=perfsuma/np.amax(perfsuma)
    plt.plot(wvl,perfsuma+cual*0.01,'b')
    plt.axvline(0,linestyle='--')
    plt.axvline(-40,linestyle='--')
    # bipos,bival=bisector(wvl,perf)
    # plt.plot(bipos,bival)
    plt.xlabel('Wavelength (km/s)')
    plt.ylabel('Intensity')
    plt.xlim(-50,25)
    plt.ylim(0.75,1.1)
    plt.subplot(1,2,2)
    
    perfobs=S['Data'][:,1]
    plt.plot(S['Data'][:,0]-27,perfobs/perfobs.max()+cual*0.01,'b')
    plt.xlim(-50,25)
    plt.ylim(0.75,1.1)
    plt.xlabel('Wavelength (km/s)')
    plt.axvline(0,linestyle='--')
    plt.axvline(-40,linestyle='--')
    #plt.ylabel('Intensity')
    plt.yticks([])
    plt.subplots_adjust(wspace=0)
    if beta==1:
        plt.savefig('Fig5_art74.png')
    if beta==.2:
        plt.savefig('Fig6_art74.png')
    
plt.show()
