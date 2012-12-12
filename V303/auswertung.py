# -*- coding: utf-8 -*-

# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.  

import math
import matplotlib

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']

def auswertung_winkel_spannung(phi, U_out, title, filename):
    (A, p), cov = opt.curve_fit(lambda x, A, p: A*np.cos(x+p), phi, U_out)

    x = np.linspace(0, 2*np.pi)

    plt.title(title)
    plt.xlabel("$\phi$")
    plt.ylabel(r"$U_\text{out}/\mathrm{V}$")
    plt.plot(phi, U_out, '+')
    plt.plot(x, A*np.cos(x+p))
    plt.grid()
    plt.savefig(filename)
    plt.close()

    print u"U_out = ({0:.5f}±{1:.5f} V)*cos(x+{2:.5f})"\
        .format(A, np.sqrt(cov[0][0]), p)

ohne_rauschen = np.loadtxt('ohne-rauschen.txt', unpack=True)
phi = np.deg2rad(ohne_rauschen[0])
U_out = ohne_rauschen[1]/1000.

auswertung_winkel_spannung(phi, U_out, "Ohne Rauschen", "ohne-rauschen.pdf")

mit_rauschen = np.loadtxt('mit-rauschen.txt', unpack=True)
phi = np.deg2rad(mit_rauschen[0])
U_out = mit_rauschen[1]/1000.

auswertung_winkel_spannung(phi, U_out, "Mit Rauschen", "mit-rauschen.pdf")

lichtmessung = np.loadtxt('messung-mit-licht.txt', unpack=True)
a = lichtmessung[0]/100 # Abstand in m
U = lichtmessung[1]/lichtmessung[2]/1000 # Spannung in V

plt.title("Die Lichtrauschmessung")
plt.xlabel("Abstand in Metern")
plt.ylabel("Spannung in Volt")
plt.grid()
plt.plot(a, U, "+")

(A, B), cov = opt.curve_fit(lambda x, A, B: A/x**2 + B, a, U)
dA = math.sqrt(cov[0][0])
dB = math.sqrt(cov[1][1])
plt.plot(a, A/a**2 + B, label=r"$\propto 1/r^2$")
plt.legend()
plt.savefig('licht-messung.pdf')

print "A = {0:.5f}+-{2:.5f}, B = {1:.5f}+-{3:.5f}".format(A, B, dA, dB)
