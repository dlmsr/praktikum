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


#1: U_K = f(I) für die Monozelle zeichnen

plot(Ig,Ug,"r*",label="Gleichspannung")




#2: Lineare Ausgleichrechnung, Steigung = -R_i, Ordinatenschnitt = U_0
(mg,yg),(dmg,dyg) = linregress.linear_fit(Ig,Ug)
(mr,yr),(dmr,dyr) = linregress.linear_fit(Ir,Ur)
(ms,ys),(dms,dys) = linregress.linear_fit(Is,Us)

def f(x,A,B):
    return A*x+B

xg= linspace(0,0.12,2)
plot(xg,f(xg,mg,yg),"r")


#Außerdem: Plot mit Gegenspannung

plot(Ic,Uc,"b*",label="Gegenspannung")
(mc,yc),(dmc,dyc) = linregress.linear_fit(Ic,Uc)
xc = linspace(0,0.12,2)
plot(xc,f(xc,mc,yc),"b")
xlabel("I / A")
ylabel("$U_K$ / V")
grid()
legend(loc="center right")
savefig("mono.pdf")
close()

# U_K = f(I) für s und r zeichnen
xr=linspace(0,0.01,2)
plot(Ir,Ur,"k*",label="Messwerte(Rechteck)")
plot(xr,f(xr,mr,yr),"k",label="Ausgleichsgerade")
xlabel("I / A")
ylabel("$U_K$ / V")
grid()
legend(loc="upper right")
savefig("rechteck.pdf")
close()

xs=linspace(0,0.002,2)
plot(Is,Us,"b*",label="Messwerte(Sinus)")
plot(xs,f(xs,ms,ys),"b",label="Ausgleichsgerade")
xlabel("I / A")
ylabel("$U_K$ / V")
grid()
legend(loc="upper right")
savefig("sinus.pdf")
close()

print("Spannungsquelle --- Innenwiderstand --- Leerlaufspannung")
print("Gleichsp. --- {} +- {} Ohm --- {} +- {} Volt".format(-mg,dmg,yg,dyg))
print("rechtecksp. --- {} +- {} Ohm --- {} +- {} Volt".format(-mr,dmr,yr,dyr))
print("Sinussp. --- {} +- {} Ohm --- {} +- {} Volt".format(-ms,dms,ys,dys))

#3: Innenwiderstand und Leerlaufsp. der Monzelle (Gegenspannung)

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
    return (U_0ges - Rges*U_0ges/(Rges + ra))**2/ra

linx = linspace(-0.1,52,1000)
plot(linx,N(linx),label = "Theoriekurve")
#Messwerte U*I gegen Ra
plot(Rg,Ug*Ig,"*", label = "Messwerte")
xlabel("$R_a$ in Ohm")
ylabel("N in Watt")
grid()
legend(loc = "center right")
savefig("leistung.pdf")
close()

print("----Leistungsbetrachtung----")
print("Theorie und Experiment stimmen überein")
