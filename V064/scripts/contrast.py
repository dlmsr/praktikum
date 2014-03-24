#! /usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import unicode_literals

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as opt

def f_max(phi, A, B):
    return A * (1 + np.sin(2 * np.deg2rad(phi) + B))

def f_min(phi, A, B):
    return A * (1 - np.sin(2 * np.deg2rad(phi) + B))

phi, I_max, I_min = np.loadtxt("../data/contrast").T

# generate figure with two axes in a row that share the x axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False)

ax2.set_xlabel("Polarisationswinkel φ im Gradmaß")
for ax in [ax1, ax2]:
    ax.set_ylabel("Spannung in Volt")

# fit data to curves
(A_max, B_max), cov = opt.curve_fit(f_max, phi, I_max)
(A_min, B_min), cov = opt.curve_fit(f_min, phi, I_min)

# generate equally distributed angles between phi_max and
# phi_min, plot data and fitted curves
alpha = np.linspace(phi.min(), phi.max())
ax1.plot(phi, I_max, "x")
ax1.plot(alpha, f_max(alpha, A_max, B_max))
ax2.plot(phi, I_min, "x")
ax2.plot(alpha, f_min(alpha, A_min, B_min))

fig.tight_layout()
fig.savefig("../figures/contrast.pdf")
