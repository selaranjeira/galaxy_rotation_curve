# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 20:51:47 2017

@author: sergio
"""

from scipy.integrate import quadrature
import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()

#======================= Definicao de cores================

OrangeLine	= "#d07000"
YellowLine	= "#d0b000"
BlueLine	= "#90a0c0"
BlueLine2       = "#483D8B"
GreyLine	= "black"
RedLine		= "#cd5c5c"
RedLine2        = "#DC143C"
PinkLine	= "#C71585"
Graycolor       = "#F5F5F5"
Seagreen        = "#20B2AA"
Turquesa        = "#40E0D0"
Roxo            = "#800080"
GreenLine       = "#008080"

#===========================================================
##---originais---###    
s = 1.32 * 5.2 * pow(10,3)       #razão massa-luminosidade do bojo: em kg/watt    
t = 3.43 * 5.2 * pow(10,3)       #razão massa-luminosidade do dissco: em Kg/watt    
s1= 1.01 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt
t1= 1.64 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt 
s2= 0.82 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt 
t2= 0.00 * 5.2 * pow(10,3)       #idem PARA PSE: em Kg/watt 

##---após o ajuste---###
ss = 3.16 * 5.2 * pow(10,3)       #razão massa-luminosidade do bojo: em kg/watt    
tt = 5.00 * 5.2 * pow(10,3)       #razão massa-luminosidade do dissco: em Kg/watt    
ss1= 1.80 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt
tt1= 0.755 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt 
ss2= 0.95 * 5.2 * pow(10,3)       #idem PARA NFW: em Kg/watt 
tt2= 2.375 * 5.2 * pow(10,3)       #idem PARA PSE: em Kg/watt                      

def v(R,rho,vel,s,t):
    
    ### ---- PARAMETROS DA GALAXIA ------ ####

    G = (6.67408 * pow(10,-11))     #unidades no SI:_m^{3}kg^{-1}s^{-2}
    
    l = 0.972 * 1.303 * pow(10,-23)     #l(0) do bojo: watt/metro^3
    k = 0.666 * 3.085 * pow(10,19)     #cte ka: metro
    N = 0.963                     #bojo: Adimensional
    q = 0.986                     #razão entre eixos: Adimensional
    
    m = 2.036 * 1.303 * pow(10,-25)   #l(0) do disco: watt/metro^3
    j = 11.304 * 3.085 * pow(10,19)     #cte ka: metro
    o = 0.564                      #N do disco: Adimensional
    u = 0.265                     #razão entre eixos - q: Adimensional
    
    #### ----- EQUACAO (5) ---- #####

    f=lambda x:(x**2)*(np.exp(-(x/k)**(1/N)))/((R**2 - (x**2)*(1-q**2))**(1./2))
    y,err = quadrature (f,0,R)

    g=lambda x:(x**2)*(np.exp(-(x/j)**(1/o)))/((R**2 - (x**2)*(1-u**2))**(1./2))
    H,err = quadrature (g,0,R)   
    
    #### ----- RESULTADO -----#####
        
    result=pow(10,-3)*(4.0*np.pi*G*(q*s*l*y + u*t*m*H + rho*vel))**(0.5)
    return (result)

####### ------CONSTRUINDO EIXO X (eixo y, result em Kms^{-1}):
erre1 = np.arange (0.05,17.5,0.05)    #em Kpc
erre = erre1 * 3.085 * pow(10,19)      #em metros

###### ---------Contribuição para vel dado perfil NFW, subindice 1 ------######
#rho1 em Kgm^{-3}
#vel1 em Kms^{-1}
r_s = 82 * 3.085 * pow(10,19)
r_s11 = 38.81 * 3.085 * pow(10,19)
rho11=(pow(10,-2.25) * 6.77 * pow(10,-20))
rho1=(pow(10,-3.72)*6.77*pow(10,-20))
vel1=(pow(r_s,3)/erre)*(np.log(1. + erre/r_s)-(erre/r_s)*(1./(1 + erre/r_s)))
vel11=(pow(r_s11,3)/erre)*(np.log(1. + erre/r_s11)-(erre/r_s11)*(1./(1 + erre/r_s11)))
########## -------parte referente ao perfil PSE, subindice 2 ----------#########
#rho2 em Kgm^{-3}
#vel2 em Kms^{-1}
r_c = 5.18 * 3.085 * pow(10,19)
r_c22 = 1.55 * 3.085 * pow(10,19)
rho22 = pow(10,-0.51)*6.77*pow(10,-20)
rho2 = pow(10,-1.24)*6.77*pow(10,-20) 
vel2 = pow(r_c,2)*(1-(r_c/erre)*np.arctan(erre/(r_c))) 
vel22 = pow(r_c22,2)*(1-(r_c22/erre)*np.arctan(erre/(r_c22))) 
########

resultados0 = np.zeros(len(erre), dtype=float)
resultados00 = np.zeros(len(erre), dtype=float)
resultados1 = np.zeros(len(erre), dtype=float)
resultados11 = np.zeros(len(erre), dtype=float)
resultados2 = np.zeros(len(erre), dtype=float)
resultados22 = np.zeros(len(erre), dtype=float)

for i in range (len(erre)):
    resultados0 [i] = v(erre[i],0,0,s,t)
    resultados00 [i] = v(erre[i],0,0,ss,tt)
    resultados1 [i] = v(erre[i],rho1,vel1[i],s1,t1)
    resultados11 [i] = v(erre[i],rho11,vel11[i],ss1,tt1)
    resultados2 [i] = v(erre[i],rho2,vel2[i],s2,t2)
    resultados22 [i] = v(erre[i],rho22,vel22[i],ss2,tt2)

####----Chi Quadrado-----####
data = np.loadtxt('/home/laranjeira/github/projects/eso509g80.dat', float)
xx, yy = data[:,0], data[:,1]
#print(len(xx))
erro1 = np.loadtxt('/home/laranjeira/github/projects/eso509g80.png.dat', float)
erroy = erro1[:,1]

erro=[]
a=int(len(erroy)/2)
for n in range(a):
    erro.append((erroy[n*2+1]-erroy[n*2])/2)


r0 = np.interp(xx,erre1,resultados0)
chi0= np.sum( np.power((yy-r0),2)/np.power(erro,2))
r00 = np.interp(xx,erre1,resultados00)
chi00= np.sum( np.power((yy-r00),2)/np.power(erro,2))
print ('azul ' ,chi0,chi00)

r1 = np.interp(xx,erre1,resultados1)
chi1= np.sum( np.power((yy-r1),2)/np.power(erro,2))
r11 = np.interp(xx,erre1,resultados11)
chi11= np.sum( np.power((yy-r11),2)/np.power(erro,2))
print ('vermelho ' ,chi1,chi11)
#plt.scatter(xx,r0,marker='*')

r2 = np.interp(xx,erre1,resultados2)
chi2= np.sum( np.power((yy-r2),2)/np.power(erro,2))
r22 = np.interp(xx,erre1,resultados22)
chi22= np.sum( np.power((yy-r22),2)/np.power(erro,2))
print ('verde ' ,chi2,chi22)

###----Plot----###
plt.plot(erre1,resultados0,c='b',label= '%.2f' % chi0)
plt.plot(erre1,resultados1,c='r',label= '%.2f' % chi1)
plt.plot(erre1,resultados2,c='g',label= '%.2f' % chi2)

plt.plot(erre1,resultados00,':b',label= '0.526')
plt.plot(erre1,resultados11,':r',label= '0.503')
plt.plot(erre1,resultados22,':g',label= '0.456')

plt.xlim(0,18.5)
plt.ylim(0,310)    
plt.legend(loc='best')
plt.errorbar(xx,yy, yerr=erro, fmt='o', markersize=1.0, color=YellowLine)
'''
#plt.text(6,60,'$r_{s}$='+str(round(((ss-s)/s)*100,2))+'% $r_{t}$='+str(round(((tt-t)/t)*100,2))+'%')
#plt.text(6,65,'$r_{s1}$='+str(round(((ss1-s1)/s1)*100,2))+'% $r_{t1}$='+str(round(((tt1-t1)/t1)*100,2))+'% $\rho_{1}$= '+ str(-3.31) + '$ \rho_{11}$= '+ str(-1.96))
#plt.text(6,70,'$r_{s2}$='+str(round(((ss2-s2)/s2)*100,2))+'% $r_{t2}$='+str(round(((tt2-t2)/t2)*100,2))+'% $\rho_{1}$= '+ str(-0.56) + '$ \rho_{11}$= '+ str(-0.77))
'''
plt.xlabel('R kpc')
plt.ylabel('v km/s')
plt.title("ESO509G80")
plt.savefig('/home/laranjeira/github/projects/ESO509G80.png')
plt.show()

elapsedtime = round(time.time()-start_time,2)
print (elapsedtime)
