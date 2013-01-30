# coding: utf-8

from __future__ import unicode_literals
import matplotlib as mpl
import numpy as np

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy import stats

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

# Planck
h = 4.135667516e-15 # eV second
# vacuum velo of light
c = 299792458 # metre per second
# diffraction distance
d = 201.4e-12 # metre


def halbwertsbreite(x, y):
    spline = UnivariateSpline(x, y-np.max(y)/2, s=0)
    r1, r2 = spline.roots() # find the roots

    lambda1 = 2*d*np.sin(np.deg2rad(r1))
    lambda2 = 2*d*np.sin(np.deg2rad(r2))
    E1 = h*c/lambda1
    E2 = h*c/lambda2
    DE = E1 - E2
    print 'Halbwertswinkel: {0:.5e} deg, {1:.5e} deg'.format(r1, r2)
    print 'Halbwertsbreite: {0:.5e}'.format(np.abs(r1-r2))
    print u'Energieaufloesung: {0:.5e} eV'.format(DE)

    xnew = np.linspace(min(x), max(x))
    ynew = spline(xnew)

    plt.plot(x, y, "x")
    plt.plot(xnew, ynew+np.max(y)/2, label='Interpolation')
    plt.axvline(r1)
    plt.axvline(r2)

    plt.grid()
    plt.legend()
    plt.xlabel("Kristallwinkel in Grad")
    plt.ylabel(u"Zählrate")

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


theta, Z = np.loadtxt("cu-emission.txt", unpack=True)

plt.subplot(211)
halbwertsbreite(theta[82:87], Z[82:87])
plt.subplot(212)
halbwertsbreite(theta[94:98], Z[94:98])
plt.savefig("halbwertsbreiten.pdf")
plt.close()

plt.xlabel("Kristallwinkel in Grad")
plt.ylabel(u"Zählrate")

plt.title("Emissionsspektrum einer Cu-Anode bei 35 kV")
plt.grid()
plt.xticks()
plt.yticks()
plt.annotate(r'$K_\alpha$', xy=(20.3, 1450))
plt.annotate(r'$K_\beta$', xy=(22.5, 5250))
plt.annotate(r'Bremsberg', xy=(10, 500))

plt.plot(theta, Z)

plt.savefig("cu-emission.pdf")
plt.close()

plot('au-absorption', "Absorptionsspektrum von Gold bei 35 kV")
plot('ge-absorption', "Absorptionsspektrum von Germanium bei 35 kV")
plot('nb-absorption', "Absorptionsspektrum von Niob bei 35 kV")
plot('rb-absorption', "Absorptionsspektrum von Rubidium bei 35 kV")
plot('zn-absorption', "Absorptionsspektrum von Zink bei 35 kV")
plot('zr-absorption', "Absorptionsspektrum von Zirkonium bei 35 kV")

theta, Z = np.loadtxt("bragg-bedingung.txt", unpack = True)
plt.plot(theta, Z)
plt.xlabel("Zählrohrwinkel in Grad")
plt.ylabel(u"Zählrate")
plt.title(u'Überprüfung der Bragg-Bedingung')
plt.grid()
plt.savefig('bragg-bedingung.pdf')
plt.close()

E = [10.9e3, 17.72e3, 14.45e3, 9.31e3, 16.89e3]
Z = [32, 41, 37, 30, 40]

slope, intercept, r_value, p_value, std_err = stats.linregress(Z, np.sqrt(E))
x = np.linspace(20, 45)

print "Moseley's law\n"
print "slope: {0:.5e}".format(slope)
print "intercept: {0:.5e}".format(intercept)
print "r-value: {0:.5e}".format(r_value)
print "p-value: {0:.5e}".format(p_value)
print "stderr: {0:.5e}".format(std_err)

plt.plot(Z, np.sqrt(E), 'x')
plt.plot(x, slope*x+intercept, 'r-', label='lin. Regression')
plt.legend()
plt.grid()
plt.title('Moseleysches Gesetz')
plt.xlabel('Kernladungszahl $Z$')
plt.ylabel('$\sqrt{E_\mathrm{K}/\mathrm{keV}}$')
plt.savefig('moseley-law.pdf')
plt.close()

