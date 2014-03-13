# encoding: utf-8

from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from scipy.constants import c #importiere Lichtgeschwindigkeit

#Messwerte laden

v_rot,C_rot,L_rot,R_rot = loadtxt("RLC_Messung_rot.txt",unpack=True)
v_schwarz,C_schwarz,L_schwarz,R_schwarz = loadtxt("RLC_Messung_schwarz.txt",unpack=True)
v_trommel,C_trommel,L_trommel,R_trommel = loadtxt("RLC_Messung_trommel.txt",unpack=True)
zeit_rot,zeit_schwarz,zeit_gruen,zeit_trommel = loadtxt("Laengenmessung.txt",unpack=True)

ampl_vorher,ampl_nachher = loadtxt("FFT.txt",unpack=True)

t1,U1 = loadtxt("k2b3.txt",unpack=True)
t2,U2 = loadtxt("k2b4.txt",unpack=True)
t3,U3 = loadtxt("k3b6.txt",unpack=True)
t1=t1*10**(-9) #s
t2=t2*10**(-9) #s
t3=t3*10**(-9) #s

v_rot = v_rot * 10**(3) #Hz
v_schwarz = v_schwarz*10**(3) #Hz
v_trommel = v_trommel*10**(3) #Hz
C_rot = C_rot * 10**(-12) #F
C_schwarz = C_schwarz * 10**(-12) #F
C_trommel = C_trommel * 10**(-9) #F
L_rot = L_rot*10**(-6) # H
L_schwarz = L_schwarz*10**(-6) #H
L_trommel = L_trommel*10**(-6) #H
zeit_rot = zeit_rot*10**(-9) #s
zeit_schwarz = zeit_schwarz*10**(-9) #s
zeit_gruen = zeit_gruen*10**(-9) #s
zeit_trommel = zeit_trommel*10**(-9) #s

#Laenge der Kabel bestimmen
def L(t):
    return c*t/(2*sqrt(2.25))

laenge_rot = L(zeit_rot) #m
laenge_schwarz = L(zeit_schwarz) #m
laenge_gruen = L(zeit_gruen) #m
laenge_trommel = L(zeit_trommel) #m

#Aus Absolutwerten die Beläge bestimmen (Absolutwert/Länge)
C_rot = C_rot/laenge_rot #F/m
C_schwarz = C_schwarz/laenge_schwarz #F/m
C_trommel = C_trommel/laenge_trommel #F/m
L_rot = L_rot/laenge_rot #H/m
L_schwarz = L_schwarz/laenge_schwarz #H/m
L_trommel = L_trommel/laenge_trommel #H/m
R_rot = R_rot/laenge_rot #Ohm/m
R_schwarz = R_schwarz/laenge_schwarz #Ohm/m
R_trommel = R_trommel/laenge_trommel #Ohm/m

#G ausrechnen
def G(R,C,L):
    return R*C/L

G_rot = G(R_rot,C_rot,L_rot) #Siemens/m
G_rot_ges = G_rot*laenge_rot #Siemens
G_schwarz = G(R_schwarz,C_schwarz,L_schwarz) #Siemens/m
G_schwarz_ges =G_schwarz*laenge_schwarz #Siemens
G_trommel = G(R_trommel,C_trommel,L_trommel) #Siemens/m
G_trommel_ges = G_trommel*laenge_trommel #Siemens

##Plot der Werte
#Widerstände
subplot(221)
plot(v_rot*10**(-3),R_rot,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),R_schwarz,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),R_trommel,'b*',label="Kabeltrommel")
semilogx()
title("Widerstandsbeläge")
xlabel("Frequenz v in kHz")
ylabel("Ohm/m")
legend(loc='center')
xlim(2.5*10**(-2),2*10*2)
grid()
#Induktivitäten
subplot(222)
plot(v_rot*10**(-3),L_rot*10**6,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),L_schwarz*10**6,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),L_trommel*10**6,'b*',label="Kabeltrommel")
semilogx()
title("Induktivitätsbeläge")
xlabel("Frequenz v in kHz")
ylabel("µH/m")
xlim(2.5*10**(-2),2*10*2)
grid()
#Kapazitäten
subplot(223)
plot(v_rot*10**(-3),C_rot*10**9,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),C_schwarz*10**9,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),C_trommel*10**9,'b*',label="Kabeltrommel")
semilogx()
title("Kapazitätsbeläge")
xlabel("Frequenz v in kHz")
ylabel("nF/m")
xlim(2.5*10**(-2),2*10*2)
grid()
#Querleitbeläge
subplot(224)
plot(v_rot*10**(-3),G_rot*10**6,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),G_schwarz*10**6,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),G_trommel*10**6,'b*',label="Kabeltrommel")
semilogx()
title("Querleitbeläge")
xlabel("Frequenz v in kHz")
ylabel("µS/m")
xlim(2.5*10**(-2),2*10*2)
grid()

tight_layout()
savefig("4er.pdf")


#Dämpfungskonstante
#Einfach FFT-Amplitude Vorher/Nachher
a_array = ampl_vorher/ampl_nachher
a = mean(a_array)
da = std(a_array)/sqrt(len(a_array))
aprom_array = a_array/laenge_trommel
aprom = mean(aprom_array)
daprom = std(aprom_array)/sqrt(len(aprom_array))

#Mehrfachreflexion
Sprung1 = 50 #V (reinkommendes Signal)
Sprung2 = 8.2 #V (an 75 Ohm reflektiert)
Sprung3 = 44.5 #V (an Ende reflektiert)

gamma_50 = Sprung2/Sprung1
gamma_75 = Sprung3/(Sprung1*(1-gamma_50))

print("Zur Mehrfachreflexion")
print("gamma_50 = {} --- gamma_75 = {}".format(gamma_50,gamma_75))

##Versch. Kästen (k) und Buchsen (b) als Abschlusswiderstände -> Messwerte fitten
#k2b3

def verlauf1(t,A,T):
    return U1[0]/2 *(2- A*(1-exp(-t/T)))    #Volt


opt1,cov1 = curve_fit(verlauf1,t1,U1)

x1=linspace(-1*10**(-8),5*10**(-7),1000)
plot(t1*10**9,U1,'r*',label="Messwerte")
plot(x1*10**9,verlauf1(x1,opt1[0],opt1[1]),'k',label="Fit")

title("Kasten 2 Buchse 3")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='upper right')
xlim(-20,475)
ylim(-3,92)
grid()
savefig("fit_k2b3.pdf")

#k2b4

def verlauf2(t,A,T):
    return U2[0]/2 *(2- A*(1-exp(-t/T)))    #Volt


opt2,cov2 = curve_fit(verlauf2,t2,U2)

x2=linspace(-1*10**(-8),5*10**(-6),1000)
plot(t2*10**9,U2,'r*',label="Messwerte")
plot(x2*10**9,verlauf2(x2,opt2[0],opt2[1]),'k',label="Fit")

title("Kasten 2 Buchse 4")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='upper right')
grid()
savefig("fit_k2b4.pdf")

#k3b6

def verlauf3(t,A,T):
    return   U3[0]/(1+A) *(2 +(A-1)*exp(-t/T))  #Volt


opt3,cov3 = curve_fit(verlauf3,t3,U3)

x3=linspace(-1*10**(-8),20*10**(-6),1000)
plot(t3*10**9,U3,'r*',label="Messwerte")
plot(x3*10**9,verlauf3(x3,opt3[0],opt3[1]),'k',label="Fit")

title("Kasten 3 Buchse 6")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='lower right')
grid()
savefig("fit_k3b6.pdf")
