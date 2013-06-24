from numpy import *
from matplotlib.pyplot import *
import linregress

####Winkelrichtgröße D bestimmen#####

(r,F,phi) = loadtxt("winkelrichtgroesse.txt",unpack = True)

r = r*10**(-2) #meter
F = F*10**(-3) #Newton
phi = deg2rad(phi) #Bogenmaß

D = 2*F*r/phi # Winkelrichtgröße in Nm

#Fehlerrechung:

deltaD = std(D)/sqrt(len(D))

print("Winkelrichtgröße= {} +- {} Nm".format(mean(D),deltaD))

####Eigenträgheitsmoment der Drillachse####
(a,T) = loadtxt("eigentraegheitsmoment.txt",unpack=True)

T = T/3 # Periodendauer in s
a = a*10**(-2) #meter
Rz=3.5/2*10**(-2) #meter
hz=3*10**(-2) #meter

m1=223.89/1000 #kg
m2=221.71/1000 #kg
m12=445.62/1000 # kg
m = (m1,m2,m12/2)
M = mean(m)
dM = mean(m)/sqrt(len(m))

#Ausgleichsgerade
(A,dA,B,dB) = linregress.linear_fit(a*a,T*T)

def f(x):
    return A*x+B

x = linspace(0,0.1,2)
grid()
plot(x,f(x),label="Ausgleichsgerade")
plot(a*a,T*T,"*",label="Messwerte")
xlabel("$a^2$ in $m^2$")
ylabel("$T^2$ in $s^2$")
legend(loc="lower right")
savefig("steiner.pdf")

Ix = mean(D)*B/(4*pi**2)-2*M*(Rz**2/4 + hz**2/12)  #-Is des Stabes!!!!


####Fehler von Ix über Gauß####

dIx=sqrt((B/(4*pi**2)*deltaD)**2 + (mean(D)/(4*pi**2)*dB)**2)
print("Eigentr. der Spirale= {} +- {} kgm^2".format(Ix,dIx))

######Trägheitsmomente von Kugel und Zylinder ######
(TK,TZ)=loadtxt("zweikoerper.txt",unpack=True)

TK=TK/4 # in s
TZ = TZ/4 # in s


#Daten der Kugel:
RK = (13.785)/2*10**(-2) #m
MK = 812.1/1000 #kg
#Daten des Zylinders:
HZ = 3*10**(-2) #m
RZ = 4.27*10**(-2)#m
MZ=1436/1000 #kg
#Theoretisch berechnete Werte:
JK_theo = 2/5*MK*RK**2
JZ_theo=MZ*RZ**2/2

def J(T):
    return T**2/(4*pi**2)*mean(D)

JK_array=J(TK)
JK=mean(JK_array)
JZ_array=J(TZ)
JZ=mean(JZ_array)

#Fehlerrechnung
dJK = sqrt(( mean(TK)*mean(D)*std(TK)/(sqrt(len(TK))*2*pi**2))**2 + (mean(TK)**2*deltaD/(4*pi**2))**2)
dJZ = sqrt(( mean(TZ)*mean(D)*std(TZ)/(sqrt(len(TZ))*2*pi**2))**2 + (mean(TZ)**2*deltaD/(4*pi**2))**2)

print("------------------------------------")
print("Trägheitsmomente zweier Körper")
print("Körper----Theoriewerte--------Messwerte")
print("Kugel-----{}---{}+-{}".format(JK_theo,JK,dJK))
print("Zylinder--{}---{}+-{}".format(JZ_theo,JZ,dJZ))
print("Verhältnis JK/JZ: Theo = {}; Mess = {}".format(JK_theo/JZ_theo,JK/JZ))



#####Puppen####
(T1,T2) = loadtxt("puppen.txt",unpack=True)
T1=T1/5#s
T2=T2/5#s

J1_array=J(T1)
J1=mean(J1_array)
J2_array=J(T2)
J2=mean(J2_array)
#Fehlerrechnung
dJ1 = sqrt(( mean(T1)*mean(D)*std(T1)/(sqrt(len(T1))*2*pi**2))**2 + (mean(T1)**2*deltaD/(4*pi**2))**2)
dJ2 = sqrt(( mean(T2)*mean(D)*std(T2)/(sqrt(len(T2))*2*pi**2))**2 + (mean(T2)**2*deltaD/(4*pi**2))**2)
#Theoriewerte
#Puppendaten
M=164.3/1000 #kg
rk=1.5*10**(-2) #m
hk=2.5*10**(-2) #m
rz=1.9*10**(-2) #m
hz=9.6*10**(-2) #m
ra=0.75*10**(-2) #m
ha=13.6*10**(-2) #m
rb=0.8*10**(-2) #m
hb=15.5*10**(-2) #m
A=2.6*10**(-2) #m

#1.Stellung:
J1_theo=M*( rk**4*hk/2 + rz**4*hz/2 + ra**4*ha + 2*ra**2*rz**2*ha + rb**4*hb + rb**2*hb*A**2/2 )/(rk**2*hk + rz**2*hz +2*ra**2*ha+2*rb**2*hb)
J2_theo=M*( rk**4*hk/2 + rz**4*hz/2 + ra**4*ha/2 + ra**2*ha**3/6 + 2*ra**2*ha*(ha/2 +rz)**2 + rb**4*hb + rb**2*hb*A**2/2 )/(rk**2*hk + rz**2*hz +2*ra**2*ha+2*rb**2*hb)

print("------------------------------------")
print("Trägheitsmomente der Puppe")
print("Stellung----Theoriewerte--------Messwerte")
print("1-----{}---{}+-{}".format(J1_theo,J1,dJK))
print("2--{}---{}+-{}".format(J2_theo,J2,dJZ))
print("Verhältnis: theo1/theo2 = {} ; prak1/prak2 = {}".format(J1_theo/J2_theo, J1/J2))
