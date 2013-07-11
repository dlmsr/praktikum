from numpy import *
from matplotlib.pyplot import *
import linregress

(Rg,Ig,Ug,Rr,Ir,Ur,Rs,Is,Us) = loadtxt("a_b_und_d.txt",unpack = True)
(Rc,Ic,Uc) = loadtxt("c.txt",unpack = True)
#Alles in SI-Einheiten umrechnen
Rg = Rg*0.5 #Ohm
Ig = Ig #Ampere
Ug = Ug # Volt

Rr = Rr*2.5 #Ohm
Ir = Ir*10**(-3)#Ampere
Ur = Ur*10**(-3)#Volt

Rs = Rs*50#Ohm
Is = Is*10**(-3)#Ampere
Us = Us*10**(-3)#Volt

Rc = Rc*0.5 #Ohm
Ic = Ic*10**(-3)#Ampere
Uc = Uc#Volt


#1: U_K = f(I), K = g,r,s zeichnen
subplot(121)
plot(Ig,Ug,"r*",label="Gleichspannung")
subplot(122)
plot(Ir,Ur,"k*",label="Rechteckspannung")
plot(Is,Us,"b*",label="Sinusspannung")


#2: Lineare Ausgleichrechnung, Steigung = -R_i, Ordinatenschnitt = U_0
(mg,yg),(dmg,dyg) = linregress.linear_fit(Ig,Ug)
(mr,yr),(dmr,dyr) = linregress.linear_fit(Ir,Ur)
(ms,ys),(dms,dys) = linregress.linear_fit(Is,Us)

def f(x,A,B):
    return A*x+B
xg=array((Ig[0],Ig[10]))
xr=array((Ir[0],Ir[10]))
xs=array((Is[0],Is[10]))
subplot(121)
plot(xg,f(xg,mg,yg),"r")
xlabel("I / A")
ylabel("$U_K$ / V")
grid()
legend(loc="upper right")
subplot(122)
plot(xr,f(xr,mr,yr),"k")
plot(xs,f(xs,ms,ys),"b")
xlabel("I / A")
ylabel("$U_K$ / V")
grid()
legend(loc="upper right")
matplotlib.pyplot.gcf().set_size_inches(16,7)
savefig("leerlauf.pdf")
close()

print("Spannungsquelle --- Innenwiderstand --- Leerlaufspannung")
print("Gleichsp. --- {} +- {} Ohm --- {} +- {} Volt".format(-mg,dmg,yg,dyg))
print("rechtecksp. --- {} +- {} Ohm --- {} +- {} Volt".format(-mr,dmr,yr,dyr))
print("Sinussp. --- {} +- {} Ohm --- {} +- {} Volt".format(-ms,dms,ys,dys))

#3: Innenwiderstand und Leerlaufsp. der Monzelle (Gegenspannung)

plot(Ic,Uc,"r*",label="Messwerte")
(mc,yc),(dmc,dyc) = linregress.linear_fit(Ic,Uc)
xc = array((Ic[0],Ic[10]))
plot(xc,f(xc,mc,yc),"r",label="Ausgleichsgerade")
grid()
xlabel("I / A")
ylabel("$U_K / V$")
legend(loc="upper left")
savefig("gegenspannung.pdf")
close()

print("--------------------------")
print("Monozelle mit Gegenspannung")
print("Innenwiederstand ------- Leerlaufspannung")
print("{} +- {} Ohm -----  {} +- {} Volt".format(mc,dmc,yc,dyc))
print("Abweichung der Mittelwerte beider Methoden")
print("Innenwiderst.={}Ohm ; Leerlaufsp.={}Volt".format(mc+mg,yg-yc))

print("Monozelle gesamt")
Rges = (-mg + mc)/2
dRges = (dmg + dmc)/2
U_0ges = (yg+yc)/2
dU_0ges = (dyg + dyc)/2
print("Innenwiderstand = {} +- {} Ohm".format(Rges,dRges))
print("Leerlaufspannung = {} +- {} Volt".format(U_0ges,dU_0ges))

#4: System.Fehler der direkten Gleichspannungsmessung

Ud = 1.5 #Volt direkt gemessen
R_v =10*10**(6) # Ohm Innenwiderstand des Messgeräts

Ud_ber = -mg # Errechnete Leerlaufspannung

Ri = Rges # Errechneter Innenwiderstand der Spannungsquelle

DUd_array = Ri*Ug/(R_v) 
DUd = mean(DUd_array)
print("-------------------------------------")
print("Systematischer Fehler der direkten Messung = {}Volt".format(DUd))

#5: Am Belastungswiderstand verrichtete Leistung der Monozelle
#Theoriekurve in Abh. von Ra

def N(ra):
    return U_0ges**2/Rges + (U_0ges - Rges*U_0ges/(Rges + ra))**2/ra

linx = linspace(Rg[0]+0.001,Rg[10])
plot(linx,N(linx),label = "Theoriekurve")
#Messwerte U*I gegen Ra
plot(Rg,Ug*Ig, label = "Messwerte")
xlabel("$R_a$ in Ohm")
ylabel("N in Watt")
grid()
legend(loc = "center right")
savefig("leistung.pdf")
close()

print("----Leistungsbetrachtung----")
print("Deutlich erkennbarer systematischer Fehler von ca. 12 Watt")
print("Es fällt deutlich weniger Leistung am Belastungswiderstand ab als errechnet")
print("Anscheinend hat das Messgerät diesen Fehler verursacht")
