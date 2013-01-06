# -*- coding: utf-8 -*-

# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.

import numpy as np
from math import pi
import matplotlib.pyplot as plt
import scipy.optimize as sco

# Daten der Apparatur
l = 632.8e-9 # Wellenlaenge des He-Ne-Lasers
L = 1 # Abstand des Spalts zum Schirm

def lade_messdaten(dateiname):
    x, I = np.loadtxt(dateiname, unpack=True)
    x = x * 1e-3 # umrechnen auf Meter
    I = I * 1e-6 # umrechnen auf A
    return (x, I)

def werte_einzelspalt_aus():
    # zu fittende Funktion Einzelspalt
    f = lambda x, A, b, c: (A*l*L * np.sin(pi*b * (x-c) / (l*L)) / (pi*(x-c)))**2

    (x, I) = lade_messdaten('single-slit.txt')
    I -= 0.8e-9 # Dunkelstrom subtrahieren

    # Fit der Meßwerte. Ersten Meßwert rausnehmen, da x=0 => div. by zero
    (A, b, c), cov = sco.curve_fit(f, x[1:], I[1:], p0=[1, 0.08e-3, 25.4e-3])

    plt.plot(x, I, "+") # Meßdaten plotten

    # Fit plotten
    x = np.linspace(0, 0.05, 1000)
    plt.plot(x, f(x, A, b, c))

    plt.title("Beugungsfigur am Einzelspalt")
    plt.xlabel("$x/\mathrm{m}$")
    plt.ylabel("$I/\mathrm{A}$")
    plt.xlim(0, 0.05)
    plt.grid()
    plt.savefig('single-slit.pdf')
    plt.close()

    print "\n Einzelspalt"
    print "A = {0:}\nb = {1:}\nc = {2:}".format(A, b, c)
    print "Fehlermatrix:\n{0:}".format(cov)


def werte_variablen_spalt_aus():
    # zu fittende Funktion Einzelspalt
    f = lambda x, A, b, c: (A*l*L * np.sin(pi*b * (x-c) / (l*L)) 
                            / (pi*(x-c)))**2

    (x, I) = lade_messdaten('variable-slit.txt')
    I -= 0.8e-9 # Dunkelstrom abziehen

    # Fit der Meßwerte. Ersten Meßwert rausnehmen, da x=0 => div. by zero
    (A, b, c), cov = sco.curve_fit(f, x[1:], I[1:],
                                   p0=[1, 0.2e-3, 25.4e-3])
    
    plt.plot(x, I, "+") # Meßdaten plotten
    
    # Fit plotten
    x = np.linspace(0, 0.05, 1000)
    plt.plot(x, f(x, A, b, c))

    plt.title("Beugungsfigur am variablen Spalt")
    plt.xlabel("$x/\mathrm{m}$")
    plt.ylabel("$I/\mathrm{A}$")
    plt.xlim(0, 0.05)
    plt.grid()
    plt.savefig('variable-slit.pdf')
    plt.close()
    
    print "\nVariabler Spalt"
    print "A = {0:}\nb = {1:}\nc = {2:}\n".format(A, b, c)
    print "Fehlermatrix:\n{0:}".format(cov)

def werte_doppelspalt_aus():
    (x, I) = lade_messdaten('double-slit.txt')
    I -= 0.8e-9 # Dunkelstrom abziehen
    b = 0.08e-3
    s = 0.33e-3
    c = 26e-3

    # zu fittende Funktion Einzelspalt
    f = lambda x, A: A * (2*l*L / (pi*b * (x-c)) 
                          * np.cos(pi*s*(x-c) / (l*L))
                          * np.sin(pi*b * (x-c) / (l*L)))**2

    # Meßdaten plotten
    plt.plot(x, I, "+") 

    # Fit plotten
    x = np.linspace(x[0]-5e-3, x[49]+5e-3, 1000)
    plt.plot(x, f(x, 0.2e-6))

    plt.title("Beugungsfigur am Doppelspalt")
    plt.xlabel("$x/\mathrm{m}$")
    plt.ylabel("$I/\mathrm{A}$")
    plt.xlim(0, 0.05)
    plt.grid()
    plt.savefig('double-slit.pdf')
    plt.close()
    
werte_einzelspalt_aus()
werte_variablen_spalt_aus()
werte_doppelspalt_aus()

