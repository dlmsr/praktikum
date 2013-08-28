# -*- coding: utf-8 -*-

# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.  

import codecs
import sys
import math
import matplotlib

from numpy import loadtxt, array_split, delete
from matplotlib import pyplot as plt
import linregress

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']


g = 9.80665 # m/s^2

def auswertung_einseitig(titel, datei, x, D, d, L, m, I, entf_x):
    u = L*x**2 - x**3/3
    entf_u = delete(u, entf_x)
    entf_D = delete(D, entf_x)
    A, dA, B, dB = linregress.linear_fit(entf_u, entf_D)

    E  = (m*g)/(2*A*I)
    dE = (m*g*dA)/(2*A**2*I)

    plt.grid()
    plt.title(titel)
    plt.xlabel(r"$(Lx^2 - \frac{x^3}{3})/\mathrm{m}^3$")
    plt.ylabel(r"$D/\mathrm{m}$")
    plt.plot(entf_u, entf_D, "r+")
    plt.plot(entf_u, A*entf_u+B)
    if len(u.take(entf_x)) != 0:
        plt.plot(u.take(entf_x), D.take(entf_x), "+", 
                 label="ausgeschlossene Werte")

        plt.legend()
    plt.savefig(datei)
    plt.clf()

    print(titel + """
  ------------------------------------------------------------------------
  Durchmesser d                     ({0:.5e} ± {1:.5e})  mm
  Elastizitätsmodul E               ({2:.5e} ± {3:.5e})  kN/mm^2
  Steigung: A                        {4:.5e} ± {5:.5e}
  ------------------------------------------------------------------------
    """.format(d.mean()*1e3, d.std(ddof=1)*1e3, E*1e-9, dE*1e-9, A, dA))

def auswertung_beidseitig(titel, datei, x, D, d, L, m, I):
    x1, x2 = array_split(x, 2)
    D1, D2 = array_split(D, 2)

    u1 = 3*L**2*x1 - 4*x1**3
    u2 = 4*x2**3 - 12*L*x2**2 + 9*L**2*x2 - L**3

    A1, dA1, B1, dB1 = linregress.linear_fit(u1, D1)
    A2, dA2, B2, dB2 = linregress.linear_fit(u2, D2)

    E1  = (m*g)/(48*A1*I)
    E2  = (m*g)/(48*A2*I)
    dE1 = (m*g*dA1)/(48*A1**2*I)
    dE2 = (m*g*dA2)/(48*A2**2*I)

    plt.title(titel)
    plt.subplot(211)
    plt.grid()
    plt.xlabel(r"$(3L^2x - 4x^3)/\mathrm{m}^3$")
    plt.ylabel(r"$D/\mathrm{m}$")
    plt.plot(u1, D1, "+")
    plt.plot(u1, A1*u1+B1)
    plt.subplot(212)
    plt.grid()
    plt.xlabel(r"$(4x^3 - 12Lx^2 + 9L^2x - L^3)/\mathrm{m}^3$")
    plt.ylabel(r"$D/\mathrm{m}$")
    plt.plot(u2, D2, "+")
    plt.plot(u2, A2*u2+B2)

    plt.savefig(datei)
    plt.clf()

    E = (E1+E2)/2

    print(titel + """
  ------------------------------------------------------------------------
  Durchmesser d                     ({0:.5e} ± {1:.5e})  mm
  Elastizitätsmodul E_f             ({2:.5e} ± {3:.5e})  kN/mm^2
  Elastizitätsmodul E_g             ({4:.5e} ± {5:.5e})  kN/mm^2
  Steigung A_f                       {6:.5e} ± {7:.5e}
  Steigung A_g                       {8:.5e} ± {9:.5e}
  Ordinatenabschnitt B_f             {10:.5e}
  Ordinatenabschnitt B_f             {11:.5e}
  ------------------------------------------------------------------------
  E-Modul (Mittelwert)               {12:.5e}  kN/mm^2
    """.format(d.mean()*1e3, d.std(ddof=1)*1e3, E1*1e-9, dE1*1e-9, E2*1e-9,
               dE2*1e-9, A1, dA1, A2, dA2, B1, B2, E))

############################################################################
# Stahl einseitig eingespannt, kreisförmiger Querschnitt

# Daten laden
x, D_v, D_n = loadtxt("stahl-einseitig-kreisfoermig.txt", unpack=True)
d = loadtxt("stahl-kreisfoermig-durchmesser.txt", unpack=True)

# Umrechnen auf Meter
x *= 1e-2
d *= 1e-3
D = (D_n - D_v)*1e-6

# Flächenträgheitsmoment
I = d.mean()**4*math.pi/64

# Auswertung beginnen
auswertung_einseitig("Stahl einseitig eingespannt, kreisförmiger Querschnitt",\
                         "stahl-einseitig-kreisfoermig.pdf",\
                         x, D, d, L = 0.495, m = 1.6809, I=I,\
                         entf_x = range(0,4)) # ersten vier Werte rausnehmen

############################################################################
# Aluminium einseitig eingespannt, quadratischer Querschnitt

# Daten laden
x, D_v, D_n = loadtxt('aluminium-einseitig-quadratisch.txt', unpack=True)
d = loadtxt('aluminium-quadratisch-durchmesser.txt', unpack=True)

# Umrechnen auf Meter
x *= 1e-2
d *= 1e-3
D = (D_n - D_v)*1e-6

# Flächenträgheitsmoment
I = d.mean()**4/12

# Auswertung beginnen
auswertung_einseitig("Aluminium einseitig eingespannt, quadratischer "\
                         "Querschnitt", "aluminium-einseitig-quadratisch.pdf",\
                         x, D, d, L = 0.496, m = 0.5416, I=I, entf_x=[])

############################################################################
# Messung beidseitig eingespannt, quadratischer Querschnitt

# Daten laden
x, D_v, D_n = loadtxt('messing-beidseitig-quadratisch.txt', unpack=True)
d = loadtxt('messing-quadratisch-durchmesser.txt', unpack=True)

# Umrechnen auf Meter
x *= 1e-2
d *= 1e-3
D = (D_n - D_v)*1e-6

# Flächenträgheitsmoment
I = d.mean()**4/12

auswertung_beidseitig("Messing beidseitig eingespannt, quadratischer "\
                          "Querschnitt", "messing-beidseitig-quadratisch.pdf",\
                          x, D, d, L=0.551, m=4.1, I=I)
