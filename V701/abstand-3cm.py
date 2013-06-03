# coding: utf-8
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from linregress import linear_fit

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
(A, B), (sA, sB) = linear_fit(x.take(values_to_use), M.take(values_to_use))

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

print("Mittlere Reichweite der Alpha-Teilchen: {0:.3e} mm".format(max))
print("Energie des Alpha-Teilchen: {0:.3e}".format((max / 3.1)**(2.0/3)))
print("Energieverlust -dE/dx = {0:.3f}+-{1:.3f}".format(-A, sA))
