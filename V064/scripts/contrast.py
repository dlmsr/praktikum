#! /usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import unicode_literals

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as opt

def f_max(phi, A, B, C):
    return A * (1 + np.sin(2 * phi + B)) + C

def f_min(phi, A, B, C):
    return A * (1 - np.sin(2 * phi + B)) + C

# printing settings
np.set_printoptions(precision=4)

phi, I_max, I_min = np.loadtxt("../data/contrast").T

# convert angles from degrees to radians
phi = np.radians(phi)

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
ax1.plot(phi, I_max, "x")
ax1.plot(alpha, f_max(alpha, popt[0,0], popt[0,1], popt[0,2]))
ax2.plot(phi, I_min, "x")
ax2.plot(alpha, f_min(alpha, popt[1,0], popt[1,1], popt[1,2]))

fig.tight_layout()
fig.savefig("../figures/contrast.pdf")

# calculate relative standard deviations of parameter estimates
std = np.zeros((2, 3))
for i in np.arange(3):
    std[:, i] = pcov[:, i, i]
std = np.sqrt(std)
relstd = std/np.abs(popt)

print("""
The estimated parameters in volts

max. intens.   A   B   C
min. intens.   A   B   C
""")
print(popt)

print("""
The relative standard deviations of the estimated parameters

max. intens.   A   B   C
min. intens.   A   B   C
""")
print(relstd)
