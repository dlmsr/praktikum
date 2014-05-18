#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#
# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.

from __future__ import unicode_literals
from __future__ import division

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize as opt
from scipy.constants import m_e, c, e, epsilon_0

# Brechungsindex von GaAs bei verwendeter Wellenlänge
n = 3.455

# Genauigkeit bei der Ausgabe
np.set_printoptions(precision=5)

z, B = np.loadtxt("../messwerte/magnetfeld", unpack=True)

plt.figure()
plt.grid()
plt.xlabel("Ort in Millimeter")
plt.ylabel("Flußdichte in Millitesla")

plt.plot(z, B, "o")

plt.savefig("../abbildungen/magnetfeld.pdf")
plt.close()

from numpy import zeros
theta1 = zeros((9,3))
theta2 = zeros((9,3))
theta  = zeros((9,3))
N      = zeros(2)
L      = zeros(2)

messwerte = np.loadtxt("../messwerte/GaAs_hochrein")

theta1[:,0] = messwerte[:,0] + messwerte[:,1]/60.0
theta2[:,0] = messwerte[:,2] + messwerte[:,3]/60.0

messwerte = np.loadtxt("../messwerte/GaAs_n1_2")

N[0] = messwerte[0,0] * 1e6 # reziproce Kubikmeter
L[0] = messwerte[1,0] * 1e-3 # Meter
theta1[:,1] = messwerte[2:,0] + messwerte[2:,1]/60.0
theta2[:,1] = messwerte[2:,2] + messwerte[2:,3]/60.0

messwerte = np.loadtxt("../messwerte/GaAs_n2_8")

N[1] = messwerte[0,0] * 1e6 # reziproce Kubikmeter
L[1] = messwerte[1,0] * 1e-3 # Meter
theta1[:,2] = messwerte[2:,0] + messwerte[2:,1]/60.0
theta2[:,2] = messwerte[2:,2] + messwerte[2:,3]/60.0

m = zeros(2)
theta   = 0.5*np.abs(theta1-theta2)
wavelen = messwerte[2:,4] # Mikrometer

print(np.column_stack((wavelen, theta1[:,0], theta2[:,0], theta[:,0],
                       theta1[:,1], theta2[:,1], theta[:,1],
                       theta1[:,2], theta2[:,2], theta[:,2])))

print(65*"-")
print("Analyse der Faraday-Rotationen\n")

fig, ax = plt.subplots(3, 1, sharex=True)
ax[2].set_xlabel("$\\lambda/\\mathrm{\mu m}$")
for i in np.arange(3):
    ax[i].plot(wavelen, theta[:,i], "o")
    ax[i].set_ylabel("$\\theta$ in Grad")
fig.tight_layout()
fig.savefig("../abbildungen/faraday-rotation.pdf")

fig, ax = plt.subplots(2, 1, sharex=True)
ax[1].set_xlabel("$\\lambda^2/\\mathrm{m}^2$")

for i in np.arange(2):
    ax[i].set_ylabel("$\\Delta\\theta$ (Bogenmaß)")

    x = (wavelen*1e-6)**2 # Quadratmeter
    y = np.radians(theta[:,i+1] - theta[:,0])

    ax[i].plot(x, y, "o")
    beta0, beta1, r_val, p_val, std_err = stats.linregress(x,y)
    print("beta0: {:.4e}\nbeta1: {:.4e}".format(beta0,beta1))
    # berechne geschätzte Standardabweichung der Parameter
    d = len(x)
    from numpy import mean, var
    s_beta0 = std_err * np.sqrt(1/d * (1 + mean(x)**2/var(x)))
    s_beta1 = std_err * np.sqrt(1/d * 1/var(x))
    # berechne effektive Masse
    m[i] = np.sqrt(e**3/8/np.pi**2/epsilon_0/c**3 * N[i]*B.max()*1e-3*L[i]/n * 1/beta0)
    print("s_beta0: {:.4e}\ns_beta1: {:.4e}".format(s_beta0,s_beta1))

    print("Effektive Masse\nm*: {:.4e}".format(m[i]))
    print("m*/m_e: {:.4f}".format(m[i]/m_e))
    x = np.linspace(x.min(), x.max())
    ax[i].plot(x, beta0*x + beta1)
    print(5*"-")

print("Mittelwerte")
print("m*: {:.4e}".format(m.mean()))
print("m*/m_e: {:.4f}".format((m/m_e).mean()))
print("Abweichung Literatur")
print("|m* - m*l|/m*l: {:.4e}".format(abs((m/m_e).mean() - 0.067)/0.067))

fig.tight_layout()
fig.savefig("../abbildungen/effektive_masse.pdf")
