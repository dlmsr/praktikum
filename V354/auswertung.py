# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from linregress import linear_fit

t, U = np.loadtxt('einhuellende.txt', unpack=True)

x = np.linspace(t.min(), t.max())
res = np.zeros((2, 4))

plt.plot(t, U, 'x')

for i in range(2):
    (m, b), err = linear_fit(t[i::2], np.log(U[i::2]))
    plt.plot(x, np.exp(m*x + b))

    res[i,0] = m
    res[i,2] = b

plt.yscale('log')
plt.grid(True, 'minor')
plt.grid(True, 'major')
plt.show()
plt.close()

print(res[:,0])

########################################################################
# Teil b)
L, C, R, R_ap = np.loadtxt('daten_und_Rap.txt', unpack=True)

R_ap_theo = np.sqrt(4*L/C)

print('Berechneter R_ap: {0:0.4e}'.format(R_ap_theo))
print('Relative Abweichung R_ap gemessen/berechnet: {0:0.4f}'\
	.format(abs(R_ap_theo - R_ap)/R_ap_theo))

########################################################################
# Teil c)
U_0 = 11 # Volt
f, U = np.loadtxt('messwerte_c.txt', unpack=True)

plt.plot(f, U/U_0, '+')
plt.grid(True, 'minor')
plt.grid(True, 'major')
plt.yscale('log')
plt.xscale('log')
plt.show()
plt.close()

########################################################################
# Teil d)
nu, deltat = np.loadtxt('messwerte_d.txt', unpack=True)

phi = 2 * np.pi * deltat * nu

plt.plot(nu, phi, '+')
plt.show()
