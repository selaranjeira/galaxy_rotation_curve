#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 28 09:44:40 2017

@author: sergio
"""

from scipy.integrate import quad
import numpy as np
#import matplotlib.pyplot as plt
import time
start_time = time.time()

def v(R,rho,vel,s,t):
    
    ### ---- PARAMETROS DA GALAXIA ------ ####

    G = (6.67408 * pow(10,-11))     #unidades no SI:_m^{3}kg^{-1}s^{-2}
    
    l = 0.972 * 1.303 * pow(10,-23)     #l(0) do bojo: watt/metro^3
    k = 0.666 * 3.085 * pow(10,19)     #cte ka: metro
    N = 0.963                     #bojo: Adimensional
    q = 0.986                     #razão entre eixos: Adimensional
    
    m = 2.036 * 1.303 * pow(10,-25)   #l(0) do disco: watt/metro^3
    j = 11.304 * 3.08 * pow(10,19)     #cte ka: metro
    o = 0.564                      #N do disco: Adimensional
    u = 0.265                     #razão entre eixos - q: Adimensional
    
    #### ----- EQUACAO (5) ---- #####

    f=lambda x:(x**2)*(np.exp(-(x/k)**(1/N)))/((R**2 - (x**2)*(1-q**2))**(1./2))
    y,err = quad (f,0,R)

    g=lambda x:(x**2)*(np.exp(-(x/j)**(1/o)))/((R**2 - (x**2)*(1-u**2))**(1./2))
    H,err = quad (g,0,R)   
    
    #### ----- RESULTADO -----#####
        
    result=pow(10,-3) * (4.0*np.pi*G*(q*s*l*y + u*t*m*H + rho*vel))**(0.5) 
    return (result)

####### ------CONSTRUINDO EIXO X (eixo y, result em Kms^{-1}):
erre1 = np.arange (0.05,18.5,0.05)          #em Kpc
erre = erre1 * 3.085 * pow(10,19)        #em metros

###### ---------Contribuição para vel dado perfil NFW, subindice 1 ------######
#rho1 em Kgm^{-3}
#vel1 em Kms^{-1}
r_s = 82 * 3.08 * pow(10,19)
#rho1=(pow(10,-5.75) * 6.77 * pow(10,-20))
vel1=(pow(r_s,3)/erre)*(np.log(1. + erre/r_s)-(erre/r_s)*(1./(1 + erre/r_s)))

########## -------parte referente ao perfil PSE, subindice 2 ----------#########
#rho2 em Kgm^{-3}
#vel2 em Kms^{-1}
r_c = 5.18 * 3.08*pow(10,19)
#rho2 = pow(10,-2.17) * 6.77 * pow(10,-20) 
vel2 = pow(r_c,2)*(1-(r_c/erre)*np.arctan(erre/(r_c))) 


#####------ERRO------#######
data = np.loadtxt('/home/laranjeira/github/projects/eso509g80.dat', float)
xx, yy = data[:,0], data[:,1]
#print(len(xx))
erro1 = np.loadtxt('/home/laranjeira/github/projects/eso509g80.png.dat', float)
erroy = erro1[:,1]

erro=[]
a=int(len(erroy)/2)


for n in range(a):
    erro.append((erroy[n*2+1]-erroy[n*2])/2)


def chinovo(y,yy):  #chi dos dados em relacao ao modelo
    chi1=np.zeros(len(y), dtype=float)
    for j in range(len(y)):
        chi1[j] = np.power((y[j]-yy[j]),2)/np.power(erro[j],2)
        #x = pow(((y[j])-yy[j])/erro[j],2)
        
    c = np.sum(chi1)
    return c

#####---------#######

#####----ajuste barionico----#####
chisquare=[]

s = 5.2 * pow(10,3)*np.linspace(2.00, 6.00, 15)        #razão massa-luminosidade do bojo: em kg/watt    
t = 5.2 * pow(10,3)*np.linspace(2.00, 6.00, 15)        #razão massa-luminosidade do dissco: em Kg/watt    

for k in range(len(s)):
    for p in range(len(t)):
        resultados0 = np.zeros(len(erre), dtype=float)
        for i in range (len(erre)):
            resultados0 [i] = v(erre[i],0,0,s[k],t[p])
        r0 = np.interp(xx,erre1,resultados0)
        xxx = chinovo(r0,yy)
        chisquare.append([s[k],t[p],xxx])


chisquare = np.array(chisquare)
res = chisquare[np.argsort(chisquare[:,2])] 
chi0 = res[0,2]
print('resultado = ', res[0,0]/(5.2 * pow(10,3)), res[0,1]/(5.2 * pow(10,3)), res[0,2])

#####-----------ajuste para NFW---------------------####
chisquare1=[]

s1 = 5.2 * pow(10,3)*np.linspace(2.00, 6.0, 15)        #razão massa-luminosidade do bojo: em kg/watt    
t1 = 5.2 * pow(10,3)*np.linspace(2.00, 6.5, 15)        #razão massa-luminosidade do dissco: em Kg/watt
alp=np.linspace(-7.0, -4.00, 15)    
rho1=pow(10,alp)* 6.77 * pow(10,-20) 
for k in range(len(s1)):
    for p in range(len(t1)):
        for a in range(len(rho1)):
            resultados1 = np.zeros(len(erre), dtype=float)
            for i in range (len(erre)):                     
                resultados1 [i] = v(erre[i],rho1[a] ,vel1[i],s1[k],t1[p])
            r1 = np.interp(xx,erre1,resultados1)
            xxx = chinovo(r1,yy)
            #print rho1 , s1[k],t1[p], xxx 
            chisquare1.append([s1[k],t1[p],rho1[a],xxx])

chisquare1 = np.array(chisquare1)
res1 = chisquare1[np.argsort(chisquare1[:,3])] 
chi1 = res1[0,3]
print('resultado1 = ', res1[0,0]/(5.2 * pow(10,3)), res1[0,1]/(5.2 * pow(10,3)), np.log10(res1[0,2]/(6.77 * pow(10,-20))), res1[0,3])


#####--------------Ajuste PSE-----######
chisquare2=[]
    
s2 = 5.2 * pow(10,3)*np.linspace(0.01, 3.0, 15)        #razão massa-luminosidade do bojo: em kg/watt    
t2 = 5.2 * pow(10,3)*np.linspace(0.01, 2.0, 15)        #razão massa-luminosidade do dissco: em Kg/watt    
bet=np.linspace(-2.00, -0.01, 10)
rho2=pow(10,bet) * 6.77 * pow(10,-20)
for k in range(len(s2)):
    for p in range(len(t2)):
        for a in range(len(rho2)):
            resultados2 = np.zeros(len(erre), dtype=float)
            for i in range (len(erre)):
                resultados2 [i] = v(erre[i],rho2[a],vel2[i],s2[k],t2[p])
            r2 = np.interp(xx,erre1,resultados2)
            xxx = chinovo(r2,yy)
            chisquare2.append([s2[k],t2[p],rho2[a],xxx])

chisquare2 = np.array(chisquare2)
res2 = chisquare2[np.argsort(chisquare2[:,3])] 
chi2 = res2[0,3]
print('resultado2 = ', res2[0,0]/(5.2 * pow(10,3)), res2[0,1]/(5.2 * pow(10,3)), np.log10(res2[0,2]/(6.77 * pow(10,-20))), res2[0,3])

####------calculo do tempo-------#####    

elapsedtime = round(time.time()-start_time,2)
print (elapsedtime)


###----Plot----###
#plt.plot(erre1,resultados0,c='b',label= '%.2f' % chi0)
#plt.plot(erre1,resultados1,c='r',label= '%.2f' % chi1)
#plt.plot(erre1,resultados2,c='g',label= '%.2f' % chi2)

#plt.plot(erre1,resultados00,':b',label= '%.2f' % chi00)
#plt.plot(erre1,resultados11,':r',label= '%.2f' % chi11)
#plt.plot(erre1,resultados22,':g',label= '%.2f' % chi22)

#plt.xlim(0,12)
#plt.ylim(0,200)    
#plt.legend(loc='best')
#plt.errorbar(xx,yy, yerr=erro , fmt='ko')
#plt.text(6,60,'$r_{s}$='+str(round(((ss-s)/s)*100,2))+'% $r_{t}$='+str(round(((tt-t)/t)*100,2))+'%')
#plt.text(6,65,'$r_{s1}$='+str(round(((ss1-s1)/s1)*100,2))+'% $r_{t1}$='+str(round(((tt1-t1)/t1)*100,2))+'%')
#plt.text(6,70,'$r_{s2}$='+str(round(((ss2-s2)/s2)*100,2))+'% $r_{t2}$='+str(round(((tt2-t2)/t2)*100,2))+'%')
#plt.xlabel('R kpc')
#plt.ylabel('v km/s')
#plt.title("ESO215G39")
#plt.show()


