# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from linregress import linear_fit

# Temp, Schwarz, Weiss, Metall rauh, Metall glänzend
data = np.loadtxt("stefan-boltzmann.txt", unpack=True)

T = data[0] + 273.16 # kelvin
U_therm = data[1:] * 1e-6 # micro volts

m = np.empty(4)
Dm = np.empty(4)
b = np.empty(4)

epsilon = np.empty(4)

for i in range(4):

    (a, c), (a_err, c_err) = linear_fit(T**4, U_therm[i])

    m[i] = a
    Dm[i] = a_err
    b[i] = c
    epsilon[i] = a/m[0]

# Gausssche Fehlerformel
Depsilon = np.sqrt((Dm/m[0])**2 + (m*Dm[0]/m[0]**2)**2)
