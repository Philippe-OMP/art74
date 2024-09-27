import numpy as np
import matplotlib.pyplot as plt
from Ylm_global import Ylm 
import pickle
from scipy.special import erf
from scipy.optimize import fsolve

def bisector(wvl,perf):

    a1=np.amin(perf)
    abajo=np.argmin(perf)
    a2=0.99#np.amax(perf)
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

cual=2
#for cual in range(30):


S=Solucion[cual]
print(S['Date'])
Y0=S['Ylm']
Star=Ylm.Create_Ylm_star(Y,Y0)
Star=Star/np.sum(abs(Star))
Starmem=Ylm.Create_Ylm_star(Y,Y0)
Vels=np.clip(Ylm.Ylm_Vels(Y,Star),-40,40)
longs,lats=np.meshgrid(Y.longitude,Y.latitude)



#plt.colorbar()
frente=(lats>0)#longs==20)
Star=Star[frente]
Vels=0.7*Vels[frente]

Acero=abs(Vels)<1
Amas=abs(Vels)>1

wvl=np.arange(-50,25,.1) 
Dopp=6.

for beta in [0.1,0.2,0.3,0.5,0.5,0.6,0.7,0.8,0.9,1]:
    perf=np.zeros(len(wvl))
    #plt.subplot(1,2,1)
    for i,w in enumerate(wvl):
        # x1=w 
        # x2=w-Vels[Amas]   
        # tau1=-(2./Vels[Amas]**2)*(np.sqrt(np.pi)/2.)*w*Dopp*(erf(x2/Dopp)-erf(x1/Dopp))    
        # tau2=-(2./Vels[Amas]**2)*(Dopp**2/2.)*(np.exp(-x2**2/Dopp**2)-np.exp(-x1**2/Dopp**2))                
        # Tau=(tau1+tau2)
        # perf[i]=np.sum(abs(Star[Amas])*(np.exp(-0.4*Tau)))+np.sum(abs(Star[Acero])*(1.-0.4*np.exp(-(w-Vels[Acero])**2/Dopp**2)))
        
        #v(z)=alfa * sqrt(1-beta z)
        #beta=.5
        x1=w-Vels[Amas]
        x2=w-Vels[Amas]*np.sqrt(1.-beta)
        tau1=(np.sqrt(np.pi)/(beta*Vels[Amas]**2))*w*Dopp*(erf(x2/Dopp)-erf(x1/Dopp))  
        tau2=(Dopp**2/(beta*Vels[Amas]**2))*(np.exp(-x2**2/Dopp**2)-np.exp(-x1**2/Dopp**2))  
        Tau=tau1+tau2 

        #v(z)=alfa (1-z)
        # x1=w-Vels[Amas]  
        # x2=w
        # tau1=np.sqrt(np.pi)*Dopp/(2.*Vels[Amas])*(erf(x2/Dopp)-erf(x1/Dopp))
        # Tau=tau1

        perf[i]=np.sum(abs(Star[Amas])*(np.exp(-0.4*Tau)))+np.sum(abs(Star[Acero])*(1.-0.4*np.exp(-(w-Vels[Acero])**2/Dopp**2)))

    perf=perf/np.amax(perf)
    #plt.plot(wvl,perf,'b')
    # plt.axvline(0,linestyle='--')
    # plt.axvline(-40,linestyle='--')
    bipos,bival=bisector(wvl,perf)
    plt.plot(bipos,bival,label=r'$\beta$ ={0}'.format(beta))
    # plt.subplot(1,2,2)
    
    # #plt.plot(S['Data'][:,0]-25,S['Data'][:,1]+cual*0.01,'b')
    # bipos,bival=bisector(S['Data'][:,0]-25,S['Data'][:,1])
    # plt.plot(bipos,bival)
    # #plt.xlim(-50,25)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Intensity')
plt.legend()   
plt.savefig('Bisectors_cual1.png')
plt.show()
