from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
import linregress

#Konstanten
c= 2.998*10**8 #m/s  Lichtgeschwindigkeit

#Messwerte laden

v_rot,C_rot,L_rot,R_rot = loadtxt("RLC_Messung_rot.txt",unpack=True)
v_schwarz,C_schwarz,L_schwarz,R_schwarz = loadtxt("RLC_Messung_schwarz.txt",unpack=True)
v_trommel,C_trommel,L_trommel,R_trommel = loadtxt("RLC_Messung_trommel.txt",unpack=True)
zeit_rot,zeit_schwarz,zeit_gruen,zeit_trommel = loadtxt("Laengenmessung.txt",unpack=True)

ampl_vorher,ampl_nachher = loadtxt("FTT.txt",unpack=True)

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
G_schwarz = G(R_schwarz,C_schwarz,L_schwarz) #Siemens/m
G_trommel = G(R_trommel,C_trommel,L_trommel) #Siemens/m

##Plot der Werte
#Widerstände
subplot(221)
plot(v_rot*10**(-3),R_rot,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),R_schwarz,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),R_trommel,'b*',label="Kabeltrommel")
title("Widerstandsbeläge")
xlabel("Frequenz v in kHz")
ylabel("Ohm/m")
legend(loc='center')
grid()
#Induktivitäten
subplot(222)
plot(v_rot*10**(-3),L_rot*10**6,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),L_schwarz*10**6,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),L_trommel*10**6,'b*',label="Kabeltrommel")
title("Induktivitätsbeläge")
xlabel("Frequenz v in kHz")
ylabel("mikroH/m")
legend(loc='center')
grid()
#Kapazitäten
subplot(223)
plot(v_rot*10**(-3),C_rot*10**9,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),C_schwarz*10**9,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),C_trommel*10**9,'b*',label="Kabeltrommel")
title("Kapazitätsbeläge")
xlabel("Frequenz v in kHz")
ylabel("nF/m")
legend(loc='center')
grid()
#Querleitbeläge
subplot(224)
plot(v_rot*10**(-3),G_rot*10**6,'r*',label="Rotes Kabel")
plot(v_schwarz*10**(-3),G_schwarz*10**6,'k*',label="schwarzes Kabel")
plot(v_trommel*10**(-3),G_trommel*10**6,'b*',label="Kabeltrommel")
title("Querleitbeläge")
xlabel("Frequenz v in kHz")
ylabel("mikroS/m")
legend(loc='center')
grid()

tight_layout()
show()


#Dämpfungskonstante
#Einfach FTT-Amplitude Vorher/Nachher
a_array = ampl_vorher/ampl_nachher
a = mean(a_array)
da = std(a_array)/sqrt(len(a_array))

#Mehrfachreflexion
Sprung1 = 50 #V (reinkommendes Signal)
Sprung2 = 15.2 #V (an 50 Ohm reflektiert)
Sprung3 = 44.5 #V (an 75 Ohm reflektiert)

gamma_50 = Sprung2/Sprung1
gamma_75 = Sprung3/Sprung1

##Versch. Kästen (k) und Buchsen (b) als Abschlusswiderstände -> Messwerte fitten
#k2b3

def verlauf1(t,T):
    return 83.25*exp(-t/T)    #Volt


opt1,cov1 = curve_fit(verlauf1,t1,U1)

x1=linspace(-1*10**(-8),5*10**(-7),1000)
plot(t1*10**9,U1,'r*',label="Messwerte")
plot(x1*10**9,verlauf1(x1,opt1[0]),'k',label="Fit")

title("Kasten 2 Buchse 3")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='upper right')
xlim(-20,475)
ylim(-3,92)
grid()
show()

#k2b4

def verlauf2(t,A,T):
    return 44.0*(A+ (1-A)*exp(-t/T))    #Volt


opt2,cov2 = curve_fit(verlauf2,t2,U2)

x2=linspace(-1*10**(-8),5*10**(-6),1000)
plot(t2*10**9,U2,'r*',label="Messwerte")
plot(x2*10**9,verlauf2(x2,opt2[0],opt2[1]),'k',label="Fit")

title("Kasten 2 Buchse 4")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='upper right')
grid()
show()

#k3b6

def verlauf3(t,A,T):
    return 10.2*(A*T+ (1-A*T)*exp(-t/T))    #Volt


opt3,cov3 = curve_fit(verlauf3,t3,U3)

x3=linspace(-1*10**(-8),20*10**(-6),1000)
plot(t3*10**9,U3,'r*',label="Messwerte")
plot(x3*10**9,verlauf3(x3,opt3[0],opt3[1]),'k',label="Fit")

title("Kasten 3 Buchse 6")
xlabel("Zeit t in ns")
ylabel("Spannung U in V")
legend(loc='upper right')
grid()
show()
