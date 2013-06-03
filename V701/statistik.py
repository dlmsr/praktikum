# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
from  matplotlib import pyplot as plt

N = np.loadtxt('statistik.txt')

print("Mittelwert der Zählraten: {0:.3f}".format(N.mean()))
print("Standardabweichung der Zählraten: {0:.3f}".format(np.std(N, ddof=1)))

# unterteile die Daten in 7 »Bins«
binnum = 7
n, low_range, binsize, extra = stats.histogram(N, binnum)

ind = np.arange(binnum)
width = 0.50

x = np.linspace(0, 7)
norm = stats.norm(4, 1.5).pdf(x)
poisson = stats.poisson(5).pmf(ind)

plt.plot(x, norm, "r", label="Normalverteilung")

plt.bar(ind, n/100., width, color="blue", label="gemessene Verteilung")

plt.bar(ind+0.5, poisson, width, color="green", label="Poisson-Verteilung")

plt.title("Statistische Auswertung des Alpha-Zerfalls")
plt.ylabel("relative Häufigkeit")

plt.xticks(ind+width, ('1', '2', '3', '4', '5', '6', '7'))
plt.grid()
plt.legend(loc="lower right")
plt.savefig("statistik.pdf")
plt.close()

