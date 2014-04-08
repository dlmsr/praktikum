#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#
# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.

from __future__ import unicode_literals
from __future__ import division

import numpy as np
from numpy import array, loadtxt, radians, linspace, zeros, mean, var
from numpy import ma
from matplotlib import pyplot as plt
from scipy import stats


# printing settings
np.set_printoptions(precision=5)

# load and prepare technical data
vac_wavelen, thickness, delta_phi, phi_1, length = loadtxt("../data/technical")

# convert lengths to nanometers
thickness = thickness * 1e6
length    = length * 1e6

# convert angles to radians
delta_phi = radians(delta_phi)
phi_1     = radians(phi_1)

########################################################################
# load data for glass measurement
counts = loadtxt("../data/refraction/glass")

# delta phi = phi_2 - phi_1
phi_2  = phi_1 + delta_phi

# calculate refraction indices and mean
ref_index = 1/(1 - counts * vac_wavelen /
               (thickness * (phi_2**2 - phi_1**2)))

print("Refraction index of glass slabs:   {:.4f} ± {:.4f}"
      .format(ref_index.mean(), ref_index.std(ddof = 1)))

########################################################################
# load and prepare data for gas measurement
pressure = zeros((11, 3))
counts   = zeros((11, 3))

pressure[:,0], counts[:,0] = loadtxt("../data/refraction/air-1", unpack=True)
pressure[:,1], counts[:,1] = loadtxt("../data/refraction/air-2", unpack=True)
pressure[:,2], counts[:,2] = loadtxt("../data/refraction/co2", unpack=True)

pressure = ma.masked_invalid(pressure)
counts   = ma.masked_invalid(counts)

# print result and plot pressure vs refraction index
ref_index = counts * vac_wavelen / length + 1

fig, ax = plt.subplots()
ax.set_xlabel("Druck in Millibar")
ax.set_ylabel("Brechungsindex")

data_points = ax.plot(pressure, ref_index, "o")
texts       = ["Luft Messung 1", "Luft Messung 2", "CO$_2$ Messung"]
plt.legend(data_points, texts, loc=4)

for i in range(3):

    # skip masked entries
    x = ma.compressed(pressure[:,i])
    y = ma.compressed(ref_index[:,i])

    beta0, beta1, r_value, p_value, std_err = stats.linregress(x, y)

    x = linspace(0, 1600)
    y = beta0 * x + beta1
    ax.plot(x, y)

    # calculate standard deviations of parameters
    n = len(x)
    s_beta0 = std_err * np.sqrt(1/n * (1 + mean(x)**2/var(x)))
    s_beta1 = std_err * np.sqrt(1/n * 1/var(x))

    print(array([[beta0, s_beta0], [beta1, s_beta1]]))
    print(array([beta0*1000 + beta1,
                 np.sqrt(1000**2*s_beta0**2 + s_beta1**2)]))

fig.savefig("../figures/pressure-ref-index.pdf")

result = np.array([pressure[:,0], ref_index[:,0], pressure[:,1],
                   ref_index[:,1], pressure[:,2], ref_index[:,2]]).T
print(result)
