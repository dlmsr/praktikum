#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

p, A, E = np.loadtxt('impulse-3cm-abstand.txt', unpack=True)

p_0 = 1013 
x = 30*p/p_0

np.where(A==24424)
max = x[16]

plt.xlabel("Effektive Weglänge in mm")
plt.ylabel("Anzahl der Impulse")
plt.title("Bestimmung der mittleren Weglänge")
plt.plot(x, A, "x")
plt.axvline(max, linestyle='--')
plt.grid()
plt.savefig('abstand-3cm.pdf')
plt.close()

values_to_use = np.arange(0, len(x)-5)

M = (max/3.1)**(2.0/3)/E[0] * E
A, B, _, _, sA = stats.linregress(x.take(values_to_use), M.take(values_to_use))

plt.title("Bestimmung des Energieverlustes $\mathrm{d}E_\\alpha/\mathrm{d}x$")
plt.xlabel("Effektive Weglänge in mm")
plt.ylabel("Energie der Teilchen in MeV")
plt.plot(x[17:], M[17:], "x", label="ausgeschlossene Werte")
plt.plot(x[:16], M[:16], "rx",
         label="Werte für die lin. Regression")
plt.plot(x[:16], A*x[:16]+B)
plt.legend()
plt.grid()
plt.savefig('abstand-3cm-energ.pdf')
plt.close()

print(f"Mittlere Reichweite der Alpha-Teilchen: {max:.3e} mm")
print(f"Energie des Alpha-Teilchen: {(max/3.1)**(2/3):.3e}")
print(f"Energieverlust -dE/dx = {-A:.3f}+-{sA:.3f}")
