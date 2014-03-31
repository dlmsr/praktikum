#! /usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import unicode_literals
from __future__ import division

import numpy as np


def calc_refraction_index_slab(thickness, vac_wavelen, counts,
                               phi_1, phi_2):

    ref_index = 1/(1 - counts * vacuum_wavelength /
                   (thickness * (phi_2**2 - phi_1**2)))

    return np.array([ref_index.mean(), ref_index.std(ddof = 1)])


def calc_refraction_index_gas(length, vac_wavelen, counts):

    ref_index = 0

    return ref_index


if __name__ == "__main__":

    # printing settings
    np.set_printoptions(precision=4)

    # load and prepare data
    data_glass        = np.loadtxt("../data/refraction/glass")
    vacuum_wavelength = data_glass[0]  # nm
    thickness         = data_glass[1]  # nm
    delta_phi         = np.radians(data_glass[2])
    phi_1             = np.radians(data_glass[3])
    counts            = data_glass[4:]

    # delta phi = phi_2 - phi_1
    phi_2  = phi_1 + delta_phi

    print(calc_refraction_index_slab(thickness, vacuum_wavelength,
                                     counts, phi_1, phi_2))

    pressure, counts = np.loadtxt("../data/refraction/air-1").T
    pressure, counts = np.loadtxt("../data/refraction/air-2").T
    pressure, counts = np.loadtxt("../data/refraction/co2").T
