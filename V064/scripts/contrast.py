#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#
# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.

from __future__ import unicode_literals

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as opt

def f_max(phi, A, B, C):
    return A * (1 + np.sin(2 * phi + B)) + C

def f_min(phi, A, B, C):
    return A * (1 - np.sin(2 * phi + B)) + C

def cov2std(cov):
    """calculate relative standard deviations of parameter estimates

    """
    std = np.zeros(len(cov))
    for i in np.arange(len(cov)):
        std[i] = cov[i, i]
    std = np.sqrt(std)
    return std


phi, I_max, I_min = np.loadtxt("../data/contrast").T

# convert angles from degrees to radians
phi = np.radians(phi)

# calculate contrast and plot data
V = (I_max - I_min)/(I_max + I_min)
x = np.linspace(0, np.pi)

fig, ax = plt.subplots()
ax.set_xlim(0, np.pi)
ax.set_xlabel("Polarisationswinkel im Bogenmaß")
ax.set_ylabel("Kontrast")
plots = ax.plot(phi, V, "o", x, np.sin(2*x))
texts = ["gemessener Kontrast", "$\sin(2x)$"]
plt.legend(plots, texts)
fig.tight_layout()
fig.show()
fig.savefig("../figures/contrast.pdf")

# generate figure with two axes in a row that share the x axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False)

ax2.set_xlim(0, np.pi)
ax2.set_xlabel("Polarisationswinkel φ im Bogenmaß")
for ax in [ax1, ax2]:
    ax.set_ylabel("Spannung in Volt")

# fit data to curves
popt = np.zeros((2, 3))      # estimated parameters
pcov = np.zeros((2, 3, 3))   # covariance matrices

popt[0,:], pcov[0,:,:] = opt.curve_fit(f_max, phi, I_max)
popt[1,:], pcov[1,:,:] = opt.curve_fit(f_min, phi, I_min)

# generate equally distributed angles between phi_max and
# phi_min, plot data and fitted curves
alpha = np.linspace(phi.min(), phi.max())
i_max = f_max(alpha, popt[0,0], popt[0,1], popt[0,2])
i_min = f_min(alpha, popt[1,0], popt[1,1], popt[1,2])

angles = np.array([phi, alpha, phi, alpha])
intens = np.array([I_max, i_max, I_min, i_min])

ax1.plot(phi, I_max, "o")
ax1.plot(alpha, i_max)
ax2.plot(phi, I_min, "o")
ax2.plot(alpha, i_min)

fig.tight_layout()
fig.savefig("../figures/intensity.pdf")
fig.show()


# printing settings
np.set_printoptions(precision=4)

print("""
The estimated parameters in volts

max. intens.   A   B   C
min. intens.   A   B   C
""")
print(popt)

std = np.zeros((2,3))
for i in np.arange(2):
    std[i,:] = cov2std(pcov[i,:,:])

relstd = std/np.abs(popt)

print("""
The relative standard deviations of the estimated parameters

max. intens.   A   B   C
min. intens.   A   B   C
""")
print(relstd)
