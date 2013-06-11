from numpy import *
from matplotlib.pyplot import *
import linregress
from scipy.optimize import curve_fit

### Guetemessung ###
(UG,fG) = loadtxt ('guetemessung.txt', unpack = True)

UE = 0.07 #mV Eingangsspannung
UG = UG/10 #Verstärkung herausgenommen

plot(fG, (UG/UE),'*', label = 'Messwerte')
grid()
xlabel('Frequenz f in kHz')
ylabel('$U_A$ / $U_E$')

#Gerade auf Höhe 1/sqrt(2) einzeichnen
x1 = [33,37]
y1 = [1/sqrt(2),1/sqrt(2)]
plot(x1,y1, label = '1/sqrt(2) - Marke')

#Gütepunkte und Filterfrequenz
plot(fG[24], (UG[24]/UE), 'r*', label = 'Gütepunkte')
plot(fG[25], (UG[25]/UE), 'r*')
plot(fG[11], (UG[11]/UE), 'r*')
plot(fG[27], (UG[27]/UE), 'r*')

filterung1 = [fG[10],fG[10]]
filterung2 = [0,(UG[10]/UE)]
plot(filterung1,filterung2,'k',label = 'Filterfrequenz')

legend(loc = 'upper left')
savefig("gueteplot.pdf")

#Berechnung der Güte Q
#Verwende dazu die Gütepunkte und die Filterfrequenz f0
f0 = fG[10]
fplus= (fG[11] + fG[27])/2
fminus = (fG[24] + fG[25])/2

Q = f0/(fplus - fminus)

print("Güte ergibt sich zu Q = {}".format(Q))
print("Keine Fehlerbestimmung möglich")

#Fehler der Güte
# Kann nicht bestimmt werden!


### Suszeptibilitätsmessung ###

(U1,R1,U2,R2,m) = loadtxt ('probenmessung.txt', unpack = True)
R1 = R1*5*10**(-3) #Ohm
R2 = R2*5*10**(-3) #Ohm
U1 = U1 * 10**(-3) /100 #Volt und Verstärkung raus
U2 = U2 * 10**(-3) /100 #Volt und Verstärkung raus

#Daten der Spule
n = 250 #Windungen
F = 86.6 * 10**(-6) # m^2 Querschnittsfläche
l = 135 * 10**(-3) # m Länge
R = 0.7 #Ohm

u0 = 4*pi*10**(-7) # N/A^2 magn. Feldkonstante
w = 35000 #Hz

#1.Methode (Neue Spannung betrachten)
#U_Sp = 1 Volt

U_1,U_2,U_3 = array_split(U2,3)

L = u0* n*n*F/l
faktor = w*w*L*L/(R*R)
print("w^2L^2 / R^2 = {}".format(faktor))
print("Verwende daher nicht die Näherungsformel, welche für w^2L^2 >> R^2 gilt")
def x1(U_Br,Q):
   return U_Br/1 *4*l/(w*u0*n*n*Q) *sqrt(R*R + w*w*(u0*n*n*F/l)**(2))

#Q_real berechnen
#Dy2O3
Q1 = 16.6/(13.5 * 7.8 ) # cm^2
Q1 = Q1 * 10**(-4) # m^2
#Nd2O3
Q2 = 9.09/(13.5 * 7.24) # cm^2
Q2 = Q2 * 10**(-4) #m^2
#Gd2O3
Q3 = 14.08/(13.5 * 7.4) # cm^2
Q3 = Q3 * 10**(-4) # m^2

susz11=x1(U_1,Q1)
susz21=x1(U_2,Q2)
susz31=x1(U_3,Q3)

#Fehlerrechnung 1. Methode
#stat.Fehler
dx11 = std(susz11)/len(susz11)
dx21 = std(susz21)/len(susz21)
dx31 = std(susz31)/len(susz31)

print("Suszeptibilitätsmessung Methode 1")
print("Dy2O3: X = {} +- {}".format(mean(susz11),dx11))
print("Nd2O3: X = {} +- {}".format(mean(susz21),dx21))
print("Gd2O3: X = {} +- {}".format(mean(susz31),dx31))
print("Fehler sind statistische Fehler")


# 2. Methode (Wiederstände neu einstellen)

def x2(dR,R3,Q):
    return 2* dR * F /((R3+998) * Q) # 998 Ohm Vorwiderstand!!!!

R_31,R_32,R_33 = array_split(R1,3)
R_neu1,R_neu2,R_neu3 = array_split(R2,3)
dR1 = R_31-R_neu1 # array mit Widerstandsdifferenzen 
dR2 = R_32-R_neu2 # Widerstand vorher - Widerstand nachher
dR3 = R_33-R_neu3 

susz12 = x2(dR1,R_31,Q1)
susz22 = x2(dR2,R_32,Q2)
susz32 = x2(dR3,R_33,Q3)

#Fehlerrechnung
#stat.Fehler
dx12 = std(susz12)/len(susz12)
dx22 = std(susz22)/len(susz22)
dx32 = std(susz32)/len(susz32)

print("Suszeptibilitätsmessung Methode 2")
print("Dy2O3: X = {} +- {}".format(mean(susz12),dx12))
print("Nd2O3: X = {} +- {}".format(mean(susz22),dx22))
print("Gd2O3: X = {} +- {}".format(mean(susz32),dx32))
