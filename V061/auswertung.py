# -*- encoding: utf-8 -*-

from __future__ import unicode_literals

from numpy import*
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

##Messwerte laden##

#Zur Wellenlänge
n1,d1 = loadtxt("messwerte/wellenlaenge/1.txt",unpack=True)
n2,d2 = loadtxt("messwerte/wellenlaenge/2.txt",unpack=True)
#Zur Polarisation
phi,poli = loadtxt("messwerte/polarisation.txt",unpack=True) # I/mA

phi = pi/180 *phi #Bogenmaß

#Zur Stabilität
L1,I1 = loadtxt("messwerte/stabilitaet/1.txt",unpack=True)
L2,I2 = loadtxt("messwerte/stabilitaet/2.txt",unpack=True)

#Zu den Moden
modpos0,modi0 = loadtxt("messwerte/mode/0.txt",unpack=True)
modpos1,modi1 = loadtxt("messwerte/mode/1.txt",unpack=True)

##Auswertung##
#Laserlicht - Wellenlänge

#Auf Abstand zum 0. Maximum umrechnen + Einheit m

d1 = d1 - 37.5
d1 = d1*10**(-2) #m
d2 = d2 - 37.6
d2 = d2*10**(-2) #m

lambda1_array = []

i=0
for i in range(0,len(n1)):
    if n1[i] != 0:
        lambda1_array.append(sin(arctan(d1[i]/0.6))/(n1[i]*100*10**(3))) #m
        i = i+1
    else:
        i = i+1

lambda2_array = []

i=0
for i in range(0,len(n2)):
    if n2[i] != 0:
        lambda2_array.append(sin(arctan(d2[i]/0.4))/(n2[i]*300*10**(3))) #m
        i = i+1
    else:
        i = i+1

i=0

lambda1= mean(lambda1_array)
lambda2 = mean(lambda2_array)
dlambda1 = std(lambda1_array)/sqrt(len(lambda1_array))
dlambda2 = std(lambda2_array)/sqrt(len(lambda2_array))
print("------Zur Wellenlänge des Lichtes-------")
print("Wellenlänge mit g=100lines/mm = {} +- {} m ".format(lambda1,dlambda1))
print("Wellenlänge mit g=300lines/mm = {} +- {} m".format(lambda2,dlambda2))
print("Wellenlänge insges. = {} +- {} m".format((lambda1 +lambda2)/2, (dlambda1+dlambda1)/2))
print()
#Polarisation des Laserlichtes

#Plot
plot(phi*180/pi,poli,'*',label="Messwerte")

a=0
for i in range(0,len(poli)):
    if poli[i] == max(poli):
        a=i
        break
    else:
        i = i+1

maxwinkel = phi[i]
i=0


print("-----Zur Polarisation des Lichtes-----")
print("Maximale Intensität bei phi = {} ° ".format(maxwinkel*180/pi))
#Besser fit durch Messwerte

def I(x,A,D):
    return A*cos(x-D)**2



fit,cov = curve_fit(I,phi,poli)
xarray=linspace(0,2*pi,100)
plot(xarray*180/pi,I(xarray,fit[0],fit[1]),'r',label="Fit")
legend(loc = 'lower right')
print("Sich durch den Fit ergebender Polarisationswinkel = {} +- {}°".format(fit[1]*180/pi,sqrt(cov[1][1])))
print()
grid()
xlabel("$\phi$ / °")
ylabel("Photostrom I / µA")
savefig("polarisation.pdf")
close()

#Stabilitätsbedingung
#(a)r1=r2=1.4m 
#(b)r1=1m ; r2=1.4m
#schöner Plot
plot(L1,I1,"-*b",label="Anordnung a")
plot(L2,I2,"-*r",label="Anordnung b")
legend(loc = 'upper left')
grid()
xlim(0,2)
xlabel("Abstand zum Laser /m")
ylabel("Photostrom I /µA")
savefig("stabilmess.pdf")
close()
print("-----Zur Stabilitätsbedingun-----")
print("a): r1=r2=1.4m :")
print("Theoretischer Stabilitätsbereich: {}m bis {}m ".format(0,2.8))
print("Ermittelter Stabilitätdbereich: {}m bis {}m".format(0,1.92))
print("b): r1=1m ;r2=1.4m :")
print("Theoretischer Stabilitätsbereich: {}m bis {} und {}m bis {}m".format(0,1.0,1.4,2.4))
print("Ermittelter Stabilitätdbereich: {}m bis {} und {}m bis {}m".format(0,1.05,1.46,1.62))
#TEM-Moden
#Grundmode

def I(r,I_0,v,w):
    return I_0*exp(-2*(r-v)**2/w**2)

modfit,modcov = curve_fit(I,modpos0,modi0)

modxarray1 = linspace(-15+modfit[1],15+modfit[1],50)
plot(modxarray1-modfit[1],I(modxarray1,modfit[0],modfit[1],modfit[2]),'-r',label="Fit")
plot(modpos0-modfit[1],modi0,'*b',label="Messwerte")
plot()
grid()
xlabel("Abs. zur opt. Achse / mm")
ylabel("Photostrom I / µA")
ylim(-0.2,1.2)
legend(loc='upper right')
savefig("mode1.pdf")
close()
#1.Oberschw.

def I2(r,I_1,I_2,v1,v2,w1,w2):
    return I_1*exp(-2*(r-v1)**2/w1**2)+ I_2*exp(-2*(r-v2)**2/w2**2)

modfit2,modcov2 = curve_fit(I2,modpos1,modi1)

minarray=linspace(5,15,1000)
minsuche = I2(minarray,modfit2[0],modfit2[1],modfit2[2],modfit2[3],modfit2[4],modfit2[5])
i=0
for i in range(0,len(minarray)):
    if minsuche[i] == min(minsuche):
        break
    else:
        i = i+1

minstelle = minarray[i]

modxarray2 = linspace(-5,27,50)
plot(modxarray2-minstelle,I2(modxarray2,modfit2[0],modfit2[1],modfit2[2],modfit2[3],modfit2[4],modfit2[5]))
plot(modpos1-minstelle,modi1,"*")
grid()
xlabel("Abs. zur opt. Achse / mm")
ylabel("Photostrom I / µA")
ylim(-2,60)
xlim(-15,15)
savefig("mode2.pdf")
close()
