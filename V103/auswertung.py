# -*- coding: utf-8 -*-

# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.  

from __future__ import unicode_literals
from __future__ import division


import codecs
import sys
import math
import matplotlib

from numpy import loadtxt, array_split, delete
from matplotlib import pyplot as plt
import linregress

UTF8Writer = codecs.getwriter('utf8')
sys.stdout = UTF8Writer(sys.stdout)

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']


g = 9.80665 # m/s^2

def auswertung_einseitig(titel, datei, x, D, d, L, m, I, entf_x):
    u = L*x**2 - x**3/3
    entf_u = delete(u, entf_x)
    entf_D = delete(D, entf_x)
    (A, B), (dA, dB) = linregress.linear_fit(entf_u, entf_D)

    E  = (m*g)/(2*A*I)
    dE = (m*g*dA)/(2*A**2*I)

    plt.grid()
    plt.title(titel)
    plt.xlabel(r"$(Lx^2 - \frac{x^3}{3})/\mathrm{m}^3$")
    plt.ylabel(r"$D/\mathrm{m}$")
    plt.plot(u.take(entf_x), D.take(entf_x), "+")
    plt.plot(entf_u, entf_D, "r+")
    plt.plot(entf_u, A*entf_u+B)
    plt.savefig(datei)
    plt.clf()

    print titel
    print 76*'-'
    print u"Durchmesser d = ({0:.3f} ± {1:.3f}) mm"\
        .format(d.mean()*1e3, d.std(ddof=1)*1e3)
    print u"Elastizitätsmodul E = ({0:.3f} ± {1:.3f}) kN/mm^2"\
        .format(E*1e-9, dE*1e-9)
    print 76*'-'

def auswertung_beidseitig(titel, datei, x, D, d, L, m, I):
    x1, x2 = array_split(x, 2)
    D1, D2 = array_split(D, 2)

    u1 = 3*L**2*x1 - 4*x1**3
    u2 = 4*x2**3 - 12*L*x2**2 + 9*L**2*x2 - L**3

    (A1, B1), (dA1, dB1) = linregress.linear_fit(u1, D1)
    (A2, B2), (dA2, dB2) = linregress.linear_fit(u2, D2)

    E1  = (m*g)/(2*A1*I)
    E2  = (m*g)/(2*A2*I)
    dE1 = (m*g*dA1)/(2*A1**2*I)
    dE2 = (m*g*dA2)/(2*A2**2*I)

    print titel
    print 76*'-'
    print u"Durchmesser d = ({0:.3f} ± {1:.3f}) mm"\
        .format(d.mean()*1e3, d.std(ddof=1)*1e3)
    print u"Elastizitätsmodul E1 = ({0:.3f} ± {1:.3f}) kN/mm^2"\
        .format(E1*1e-9, dE1*1e-9)
    print u"Elastizitätsmodul E2 = ({0:.3f} ± {1:.3f}) kN/mm^2"\
        .format(E2*1e-9, dE2*1e-9)
    print 76*'-'

############################################################################
# Stahl einseitig eingespannt, kreisförmiger Querschnitt

# Daten laden
x, D_v, D_n = loadtxt('stahl-einseitig-kreisfoermig.txt', unpack=True)
d = loadtxt('stahl-kreisfoermig-durchmesser.txt', unpack=True)

# Umrechnen auf Meter
x *= 1e-2
d *= 1e-3
D = (D_n - D_v)*1e-6

# Flächenträgheitsmoment
I = d.mean()**4*math.pi/64

# Auswertung beginnen
auswertung_einseitig(u"Stahl einseitig eingespannt, kreisförmiger Querschnitt",\
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
auswertung_einseitig(u"Aluminium einseitig eingespannt, quadratischer "\
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

auswertung_beidseitig(u"Messing beidseitig eingespannt, quadratischer "\
                          "Querschnitt", "messing-beidseitig-quadratisch.txt",\
                          x, D, d, L=55.1, m=4.1, I=I)
