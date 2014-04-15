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

#Polarisation des Laserlichtes

#Plot
plot(phi,poli,'*',label="Messwerte")

a=0
for i in range(0,len(poli)):
    if poli[i] == max(poli):
        a=i
        break
    else:
        i = i+1

maxwinkel = phi[i]
i=0


print("Zur Polarisation des Lichtes")
print("Maximale Intensität bei phi = {} rad".format(maxwinkel))
#Besser fit durch Messwerte

def I(x,A,D):
    return A*cos(x-D)**2



fit,cov = curve_fit(I,phi,poli)
xarray=linspace(0,2*pi,100)
plot(xarray,I(xarray,fit[0],fit[1]),'r',label="Fit")
legend(loc = 'lower right')
print("Sich durch den Fit ergebender Polarisationswinkel = {} rad".format(fit[1]))
grid()
show()

#

#Stabilitätsbedingung
#r1=r2=1.4m

#schöner Plot
plot(L1,I1,"b",label="1")
plot(L2,I2,"r",label="2")
legend(loc = 'upper right')
grid()
show()

#TEM-Moden
#Grundmode
plot(modpos0,modi0,"*")
grid()
show()
#1.Oberschw.
plot(modpos1,modi1,"*")
grid()
show()
