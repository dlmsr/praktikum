# encoding: utf-8
from numpy import *
from matplotlib.pyplot import *
import linregress

# ########### Teil 1 (Linse 1 mit f = 50 mm) und 2 (Wasserlinse) ############

(g_1, b_1 , g_2 , b_2) = loadtxt('messwerte_linsen1und2.txt', unpack = True) # alles in cm

#Brennweiten bestimmen.

def f(g,b):
    return (g*b)/(g+b)


f_1 = f(g_1,b_1) # cm
f_2 = f(g_2,b_2) # cm

f_1 = f_1 * 10 # in mm umrechnen
f_2 = f_2 * 10 # in mm umrechnen

# Fehlerrechnung:

df_1 = std(f_1)/len(f_1) # mm
df_2 = std(f_2)/len(f_2) # mm

# Überprüfung der Güte der Messung anhand eines Plots

#Nullarray zum Plotten der g's und b's auf den x.- bzw. y-Achsen:
nullen = zeros(len(g_1))

plot(g_1,nullen, "*")
plot(nullen,b_1,"*") 

#Geraden durch die Punkte
def g(x,m,b):
    return m*x +b

m_1 = -b_1 / g_1

i = 1

# Ein Geradenplot wird außerhalb der Schleife durchgeführt, um das Label für die Legende zu setzen!
I = linspace(0 , g_1[1], 2)
plot(I, g(I,m_1[1],b_1[1]), 'b', label = 'Linse 1 ( f = 5 cm)')

for i in range(1 , len(g_1)):
   I = linspace(0 , g_1[i], 2)
   plot(I, g(I,m_1[i],b_1[i]), 'b') 


plot([4.84], [4.84], 'g.', markersize=15.0)

# Das gleiche mit der Wasserlinse
nullen = zeros(len(g_1))

plot(g_2,nullen, "*")
plot(nullen,b_2,"*")

#Geraden durch die Punkte
m_2 = -b_2 / g_2

i = 1

# Ein Geradenplot wird außerhalb der Schleife durchgeführt, um das Label für die Legende zu setzen!
I = linspace(0 , g_2[0], 2)
plot(I, g(I,m_2[0],b_2[0]), 'r', label = 'Linse 2 ( f unbekannt)')
   
for i in range(1 , len(g_2)):
   I = linspace(0 , g_2[i], 2)
   plot(I, g(I,m_2[i],b_2[i]), 'r') 


plot([10.283], [10.283], 'g.', markersize=15.0, label = 'Berechnete Brennweiten')


legend(loc = 'upper right')
xlabel('Gegenstandsweite g/cm')
ylabel('Bildweite b/cm')
xlim([0,45])
ylim([0,45])
grid()
savefig("linsen1und2.pdf")
close()

##################### Teil 3 & 4 #############################
# Bessel für weißes, rotes und blaues Licht

(e_1,b2w,g2w,g1w,b1w,e_20,b2r0,g2r0,g1r0,b1r0,b2b0,g2b0,g1b0,b1b0) = loadtxt('messwerte_bessel.txt',unpack = True)

# arrays für rot und blau sind zur hälfte mit Nullen gefüllt. Arrays richtig machen:
(e_2 , e_2muell) = array_split(e_20, 2)
(b2r , b2rmuell) = array_split(b2r0 ,2)
(g2r , g2rmuell) = array_split(g2r0, 2)
(g1r , g1rmuell) = array_split(g1r0, 2)
(b1r , b1rmuell) = array_split(b1r0, 2)
(b2b , b2bmuell) = array_split(b2b0 ,2)
(g2b , g2bmuell) = array_split(g2b0, 2)
(g1b , g1bmuell) = array_split(g1b0, 2)
(b1b , b1bmuell) = array_split(b1b0, 2)


def fb(e,d):
    return (e**2 - d**2)/(4*e)

# 1.: weißes Licht
ew = e_1
dw1 = abs(g1w - b1w)
dw2 = abs(g2w - b2w)
fw1 = fb(ew,dw1) #cm
fw2 = fb(ew, dw2)

# 2.: rotes Licht
er = e_2
dr1 = abs(g1r - b1r)
dr2 = abs(g2r - b2r)
fr1 = fb(er, dr1)
fr2 = fb(er, dr2)

# 3.: blaues Licht
eb = e_2
db1 = abs(g1b - b1b)
db2 = abs(g2b - b2b)
fb1 = fb(eb,db1)
fb2 = fb(eb,db2)

#Fehlerrechnung

deltafw1 = std(fw1)/len(fw1)
deltafw2 = std(fw2)/len(fw2)

deltafr1 = std(fr1)/len(fr1)
deltafr2 = std(fr2)/len(fr2)

deltafb1 = std(fb1)/len(fb1)
deltafb2 = std(fb2)/len(fb2)

###################### Teil 5 ############################
# Abbe: Linsensystem mit f= -100 und f = 100
(gs,bs,B) = loadtxt('messwerte_abbe.txt', unpack = True) # in cm
G = 3.8 # cm
V = B/G

plot((1+ 1/V), gs, 'b*', label = 'Zu den Gegenstandsweiten')

#Ausgleichsgerade
def gerade(x,m,y):
    return m*x+y

(M1,Y1),cov1 = linregress.linear_fit((1+ 1/V), gs)
u = linspace(min(V),max(V),2)
plot(1+ 1/u, gerade(1 + 1/u,M1,Y1), "b")


#bs
plot((1+ V), bs,"r*", label = 'Zu den Bildweiten')
(M2,Y2),cov2 = linregress.linear_fit((1+V), bs)
u = linspace(min(V),max(V),2)
plot(1+u, gerade(1 + u,M2,Y2), "r")
ylabel(' cm')
grid()
legend(loc = 'lower right')
savefig("abbe.pdf")
close()

#Fehlerrechnung

deltafg = cov1[0] #cm
deltafb = cov2[0] #cm

deltah = cov1[1] #cm
deltahs = cov2[1] #cm
