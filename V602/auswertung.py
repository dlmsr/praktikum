# coding: utf-8

from __future__ import unicode_literals
import matplotlib as mpl
import numpy as np

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

theta, Z = np.loadtxt("cu-emission.txt", unpack=True)

x = theta[82:87]
y = Z[82:87]

spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
r1, r2 = spline.roots() # find the roots

print r1, r2

xnew = np.linspace(min(x), max(x))
ynew = spline(xnew)

plt.plot(x, y, "x")
plt.plot(xnew, ynew+np.max(y)/2)

plt.xlabel("Kristallwinkel in Grad")
plt.ylabel(u"Zählrate")
plt.title(r"Halbwertsbreite der $K_\alpha$-Linie der Cu-Anode bei 35 kV")
plt.grid(which='both')
plt.show()
plt.close()

def plot(dateiname, titel):
    theta, Z = np.loadtxt(dateiname+".txt", unpack=True)

    plt.xlabel("Kristallwinkel in Grad")
    plt.ylabel(u"Zählrate")

    plt.title(titel)
    plt.grid()
    plt.xticks()
    plt.yticks()

    plt.plot(theta, Z)

    plt.savefig(dateiname+".pdf")
    plt.close()

plot('cu-emission', "Emissionsspektrum einer Cu-Anode bei 35 kV")
plot('au-absorption', "Absorptionsspektrum von Gold bei 35 kV")
plot('ge-absorption', "Absorptionsspektrum von Germanium bei 35 kV")
plot('nb-absorption', "Absorptionsspektrum von Niob bei 35 kV")
plot('rb-absorption', "Absorptionsspektrum von Rubidium bei 35 kV")
plot('zn-absorption', "Absorptionsspektrum von Zink bei 35 kV")
plot('zr-absorption', "Absorptionsspektrum von Zirkonium bei 35 kV")

