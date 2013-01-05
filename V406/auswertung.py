# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License. To view a copy of this
# license, visit http://creativecommons.org/licenses/by-nc/3.0/.

import numpy as np
from math import pi
import matplotlib.pyplot as plt
import scipy.optimize as sco

def lade_messdaten(dateiname):
    x, I = np.loadtxt(dateiname, unpack=True)
    x *= 1e-3 # umrechnen auf Meter
    I *= 1e-6 # umrechnen auf A
    return (x, I)

l = 632.8e-9 # Wellenlaenge des He-Ne-Lasers
L = 1 # Abstand des Spalts zum Schirm

# zu fittende Funktion Einzelspalt
f = lambda x, A, b, c: (A * l * L * np.sin(pi*b * (x-c) / l / L) / pi / (x-c))**2

(x, I) = lade_messdaten('single-slit.txt')
I -= 0.8e-9 # Dunkelstrom subtrahieren

# Fit der Meßwerte. Ersten Meßwert rausnehmen, da x=0 => div. by zero
(A, b, c), cov = sco.curve_fit(f, x[1:], I[1:], p0=[1, 0.08e-3, 25.4e-3])

plt.plot(x, I, "+") # Meßdaten plotten

# Fit plotten
x = np.linspace(x[0]-5e-3, x[49]+5e-3, 1000)
plt.plot(x, f(x, A, b, c))

print A
print b
print c
print cov

plt.grid()
plt.show()

(x, I) = lade_messdaten('variable-slit.txt')
I -= 0.8e-9 # Dunkelstrom abziehen

# Fit der Meßwerte. Ersten Meßwert rausnehmen, da x=0 => div. by zero
(A, b, c), cov = sco.curve_fit(f, x[1:], I[1:], p0=[1, 0.2e-3, 25.4e-3])

plt.plot(x, I, "+") # Meßdaten plotten

# Fit plotten
x = np.linspace(x[0]-5e-3, x[49]+5e-3, 1000)
plt.plot(x, f(x, A, b, c))

print A
print b
print c
print cov

plt.grid()
plt.show()

(x, I) = lade_messdaten('double-slit.txt')
I -= 0.8e-9 # Dunkelstrom abziehen

# zu fittende Funktion Einzelspalt
f = lambda x, A, b, c: 4*cos(pi)

# Fit der Meßwerte. Ersten Meßwert rausnehmen, da x=0 => div. by zero
(A, b, c), cov = sco.curve_fit(f, x[1:], I[1:], p0=[1, 0.2e-3, 25.4e-3])

plt.plot(x, I, "+") # Meßdaten plotten

# Fit plotten
x = np.linspace(x[0]-5e-3, x[49]+5e-3, 1000)
plt.plot(x, f(x, A, b, c))

print A
print b
print c
print cov

plt.grid()
plt.show()
