# -*- coding: utf-8 -*-

# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.  

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

def auswertung_winkel_spannung(phi, U_out, p0=None):
    (A, p), cov = opt.curve_fit(lambda x, A, p: A*np.cos(x+p), phi, U_out, p0)

    x = np.linspace(0, 2*np.pi)
    plt.plot(phi, U_out, '+')
    plt.plot(x, A*np.cos(x+p))
    plt.show()

    print u"U_out = ({0:.5f}±{1:.5f} V)*cos(x+{2:.5f})"\
        .format(A, np.sqrt(cov[0][0]), p)

ohne_rauschen = np.loadtxt('ohne-rauschen.txt', unpack=True)
phi = np.deg2rad(ohne_rauschen[0])
U_out = ohne_rauschen[1]/1000.

auswertung_winkel_spannung(phi, U_out)

mit_rauschen = np.loadtxt('mit-rauschen.txt', unpack=True)
phi = np.deg2rad(mit_rauschen[0])
U_out = mit_rauschen[1]/1000.

auswertung_winkel_spannung(phi, U_out)
