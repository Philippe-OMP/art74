import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf,roots_legendre
from scipy.optimize import fsolve
from datetime import datetime

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
plt.subplot(1,2,1)
core=[]
cont=[]
fechas=[]
for cual in range(cuantos-6):

    S=Solucion[cual]
    print(S['Date'])
    fechas.append(datetime.strptime(S['Date'],'%Y/%m/%d'))
    Y0=S['Ylm']
    Star=Ylm.Create_Ylm_star(Y,Y0)
    Star=Star/np.sum(abs(Star))
    Starmem=Ylm.Create_Ylm_star(Y,Y0)
    Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-40,10)
    longs,lats=np.meshgrid(Y.longitude,Y.latitude)
    
    wvl=np.arange(-50,25,.1) 
    nodos=30
    roots,weights=roots_legendre(nodos) 
    phis=np.pi*np.arange(nodos+1)/nodos-np.pi/2
    thetas=np.arccos(roots)     
    beta=1
    Dopp=6
    perfsuma=np.zeros(len(wvl))
    photo=0
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
            photo=photo+abs(Star[estelat,estelong])

    print('Fotometria',photo)
    
    cont.append(np.amax(perfsuma))
    # plt.plot(wvl,perfsuma)
    # plt.show()
    perfsuma=perfsuma/np.amax(perfsuma)
    core.append(wvl[np.argmin(perfsuma)])

cont=np.asarray(cont)
cont=cont/cont.max()

plt.plot(cont,core,'-o')
plt.xlim(0.7,1.03)
plt.xlabel('Continuum brightness')
plt.ylabel('Line core velocity (km/s)')
plt.subplot(1,2,2)
plt.plot(fechas,core,'o')
plt.xlabel('Date')
plt.yticks([])
plt.subplots_adjust(wspace=0)
plt.savefig('Fig11_art74.png')

plt.show()
