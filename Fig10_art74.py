import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf,roots_legendre
from scipy.optimize import fsolve

def bisector(wvl,perf):

    a1=np.amin(perf)
    abajo=np.argmin(perf)
    a2=0.98#np.amax(perf)
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






def RT(wvl,vcos,Dopp,beta):
    x= lambda z: np.sqrt(1.-beta*z)
 
    if abs(vcos)>1:               
        tau1=np.sqrt(np.pi)*wvl*Dopp*(erf((wvl-vcos*x(1))/Dopp)-erf((wvl-vcos*x(0))/Dopp))  
        tau2=Dopp**2*(np.exp(-(wvl-vcos*x(1))**2/Dopp**2)-np.exp(-(wvl-vcos*x(0))**2/Dopp**2))  
        Tau=(tau1+tau2)/(beta*vcos**2)   # This vcos**2 is what justifies the separation of cases for vcos ~ 0
        
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    perfil=np.exp(-0.7*Tau)
    return perfil

def RT_linear(wvl,vcos,Dopp,beta):
  
    if abs(vcos)>1:
        Tau=-np.sqrt(np.pi)*Dopp*(erf((wvl-vcos*beta)/Dopp)-erf((wvl/Dopp)))    
        Tau=Tau/(beta*vcos*2)
    else:
        Tau=np.exp(-(wvl-vcos)**2/Dopp**2)
    Tau=-np.sqrt(np.pi)*Dopp*(erf((wvl-vcos*beta)/Dopp)-erf((wvl/Dopp)))    
    Tau=Tau/(beta*vcos*2)
    # print('Tau:',Tau.max(),Tau.min())
    # if Tau.min()<-0.01:
    #     print(vcos,Dopp,beta)
    perfil=np.exp(-Tau)
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

cual=2
#for cual in range(30):


S=Solucion[cual]
Y0=S['Ylm']
Star=Ylm.Create_Ylm_star(Y,Y0)
Star=Star/np.sum(abs(Star))
Starmem=Ylm.Create_Ylm_star(Y,Y0)
Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-40,40)
longs,lats=np.meshgrid(Y.longitude,Y.latitude)

wvl=np.arange(-50,25,.1) 
nodos=20
roots,weights=roots_legendre(nodos) 
phis=np.pi*np.arange(nodos+1)/nodos-np.pi/2
thetas=np.arccos(roots)     
perf=np.zeros(len(wvl))
Dopp=6

for beta in [0.1,0.5,1,1.5,2.]:
    perfsuma=np.zeros(len(wvl))
    for phi in phis:
        for i,theta in enumerate(thetas):
            lat=coord2lat(theta,phi)
            estelat=np.argmin(abs(Y.latitude*rad-lat))
            az=coord2long(theta,phi)
            estelong=np.argmin(abs(Y.longitude*rad-az))
            #print(abs(Star[estelat,estelong]),Vels[estelat,estelong],np.sin(lat))
            #print(phi/rad,theta/rad,lat/rad,az/rad)
            perf=abs(Star[estelat,estelong])*RT_linear(wvl,Vels[estelat,estelong],Dopp,beta)*np.sin(theta)*np.cos(phi)
            perfsuma=perfsuma+(np.pi/nodos)*weights[i]*perf

    perfsuma=perfsuma/np.amax(perfsuma)
    
    # plt.plot(wvl,perfsuma,'b')
    # plt.axvline(0,linestyle='--')
    # plt.axvline(-40,linestyle='--')
    bipos,bival=bisector(wvl,perfsuma)
    plt.plot(bipos,bival,label=r'$\beta$ ={0}'.format(beta))
    # plt.subplot(1,2,2)
    
    # #plt.plot(S['Data'][:,0]-25,S['Data'][:,1]+cual*0.01,'b')
    # bipos,bival=bisector(S['Data'][:,0]-25,S['Data'][:,1])
    # plt.plot(bipos,bival)
    # #plt.xlim(-50,25)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Intensity')
plt.legend()  
plt.savefig('Fig10_art74.png') 
plt.show()
