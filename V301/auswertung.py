import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

plt.style.use("mypubstyle")
plt.rcParams.update(**{
    "pgf.preamble": [r"\usepackage{siunitx}"],
    "text.latex.preamble": [r"\usepackage{siunitx}"],
})

(Rg,Ig,Ug,Rr,Ir,Ur,Rs,Is,Us) = np.loadtxt("a_b_und_d.txt",unpack = True)
(Rc,Ic,Uc) = np.loadtxt("c.txt",unpack = True)
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

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(5.78, 2.2),
                                    squeeze=True)

#1: U_K = f(I) für die Monozelle zeichnen

ax1.plot(Ig,Ug,"r*",label="Gleichspannung")




#2: Lineare Ausgleichrechnung, Steigung = -R_i, Ordinatenschnitt = U_0
mg, yg, _, _, dmg = stats.linregress(Ig, Ug)
mr, yr, _, _, dmr = stats.linregress(Ir, Ur)
ms, ys, _, _, dms = stats.linregress(Is, Us)

dyg = dmg * np.sqrt(1/(Ig.size) * np.sum(Ig**2))
dyr = dmr * np.sqrt(1/(Ir.size) * np.sum(Ir**2))
dys = dms * np.sqrt(1/(Is.size) * np.sum(Is**2))


def f(x,A,B):
    return A*x+B

xg= np.linspace(0,0.12,2)
ax1.plot(xg,f(xg,mg,yg),"r")


#Außerdem: Plot mit Gegenspannung

ax1.plot(Ic,Uc,"b*",label="Gegenspannung")
mc, yc, _, _, dmc = stats.linregress(Ic, Uc)
dyc = dmc * np.sqrt(1/(Ic.size) * np.sum(Ic**2))

xc = np.linspace(0,0.12,2)
ax1.plot(xc,f(xc,mc,yc),"b")
ax1.set_xlabel(r"$I/\si{\ampere}$")
ax1.set_ylabel(r"$U_K/\si{\volt}$")
ax1.set_title("Monozelle")
ax1.legend(loc="upper left", fontsize=6)


# U_K = f(I) für s und r zeichnen
ax2.set_title("Rechteck")
ax2.set_xlabel(r"$I/\si{\ampere}$")
xr=np.linspace(0,0.01,2)
ax2.plot(Ir,Ur,"k*")
ax2.plot(xr,f(xr,mr,yr),"k")


ax3.set_title("Sinus")
ax3.set_xlabel(r"$I/\si{\ampere}$")
xs=np.linspace(0,0.002,2)
ax3.plot(Is,Us,"k*")
ax3.plot(xs,f(xs,ms,ys),"k")


plt.tight_layout(pad=0.25)
fig.savefig("linregress.pdf")

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
DUd = np.mean(DUd_array)
print("-------------------------------------")
print("Systematischer Fehler der direkten Messung = {}Volt".format(DUd))

#5: Am Belastungswiderstand verrichtete Leistung der Monozelle
#Theoriekurve in Abh. von Ra

def N(ra):
    return (U_0ges - Rges*U_0ges/(Rges + ra))**2/ra

linx = np.linspace(-0.1,52,1000)

plt.figure(figsize=(5.78, 2.89))
plt.plot(linx, N(linx), color="black", label="Theoriekurve")
#Messwerte U*I gegen Ra
plt.plot(Rg,Ug*Ig, "k*", label="Messwerte")
plt.xlabel("$R_a$ in Ohm")
plt.ylabel("N in Watt")
plt.grid()
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("leistung.pdf")
plt.close()

print("----Leistungsbetrachtung----")
print("Theorie und Experiment stimmen überein")
