# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from linregress import linear_fit

# Laden der Geräte-Daten
L, C, R, R_ap = np.loadtxt('daten_und_Rap.txt', unpack=True)

########################################################################
# Teil a)

# Lade Daten und lege Ergebnis-Matrix an
minmax         = np.loadtxt('einhuellende.txt')
minima, maxima = np.vsplit(minmax.T, 2)
res            = []

for (t, U) in [maxima, minima]:
    # Plotte Meßwerte
    plt.semilogy(t, U, 'x')

    # Berechne Regressionsgerade und speichere Werte
    m, sm, b, sb = linear_fit(t, np.log(U))
    res.append([m, sm, b, sb])

    # Plotte Regressionsgerade
    x = np.linspace(t.min(), t.max())
    plt.semilogy(x, np.exp(m*x + b))

# konvertiere in numpy array
res = np.array(res)

plt.title('Bestimmung des Exponenten')
plt.ylabel('$U_C/\\mathrm{V}$')
plt.xlabel('$t/\\mathrm{s}$')
plt.ticklabel_format(axis='x',style='sci',scilimits=(1,4))
plt.grid(True, which='both')

plt.savefig('exp-plot.pdf')
plt.close()

# Mittelwert und Standardabweichung vom Mittelwert der Steigungen
res = np.append(res, [np.mean(res, axis = 0), stats.sem(res, axis = 0)],
                axis = 0)

# effektiver Dämpfungswiderstand
R_eff = -2*res[2,0]*L

print('Ergebnis der Ausgleichsrechnung:')
print('m\tb\tsm\tsb')
print(res)
print('{: >8.4f}\teffektiver Dämpfungswiderstand'.format(R_eff))
print('{: >8.4f}\teingebauter Widerstand'.format(R))
print('{: >8.2%}\tRelative Abweichung vom berechneten Wert'.format(abs(R - R_eff)/R))
print(72*'-')

########################################################################
# Teil b)

# Berechne R_ap aus Gerätedaten
R_ap_theo = 2*np.sqrt(L/C)

print('Berechneter R_ap: \t{:0.4e}'.format(R_ap_theo))
print('Relative Abweichung: \t{:0.4%}'\
	.format(abs(R_ap_theo - R_ap)/R_ap_theo))
print(72*'-')

########################################################################
# Teil c)
U_0 = 11 # Volt
f, U = np.loadtxt('messwerte_c.txt', unpack=True)

plt.title('Bestimmung der Resonanzüberhöhung')

plt.subplot(211)
plt.semilogy(f, U/U_0, '+')
plt.ylabel('$U/U_0$')
plt.xlabel('$\\nu/\\mathrm{Hz}$')
plt.grid(True, which='both')
plt.ticklabel_format(axis='x',style='sci',scilimits=(1,4))

sec = np.where(U/U_0 >= 5)
plt.subplot(212)
plt.plot(f[sec], U[sec]/U_0, '+')
plt.ylabel('$U/U_0$')
plt.xlabel('$\\nu/\\mathrm{Hz}$')
plt.grid(True, which='both')
plt.ticklabel_format(axis='x',style='sci',scilimits=(1,4))

plt.savefig('resonanz.pdf')
plt.close()

q_theo = 1/R*np.sqrt(L/C)
q = U.max()/U_0
breite_theo = R/L
breite = (f[sec].max() - f[sec].min())

print('Berechnete Resonanzüberhöhung: {0:0.4e}'.format(q_theo))
print('Experimentell bestimmte Resonanzüberhöhung: {0:0.4e}'.format(q))
print('relative Abweichung: {0:0.4%}'.format(abs(q_theo-q)/q_theo))
print('Berechnete Breite der Resonanzkurve: {0:0.4e}'.format(breite_theo))
print('Experimentell bestimmte Breite: {0:0.4e}'.format(breite))
print('relative Abweichung: {0:0.4%}'\
      .format(abs(breite_theo-breite)/breite_theo))
print(72*'-')

########################################################################
# Teil d)
nu, deltat = np.loadtxt('messwerte_d.txt', unpack=True)

phi = 2 * np.pi * deltat * 1e-9 * nu * 1e3

plt.subplot(211)
plt.semilogy(nu, phi, '+')
plt.grid(True, which='both')
plt.xlabel('$\\nu/\\mathrm{kHz}$')
plt.ylabel('$\\phi$')

sec = np.where(nu > 25)

plt.subplot(212)
plt.grid(True, which='both')
plt.xlabel('$\\nu/\\mathrm{kHz}$')
plt.ylabel('$\\phi$')
plt.plot(nu[sec], phi[sec], '+')

plt.title('Phase gegen Frequenz')
plt.savefig('phasen-plot.pdf')
plt.close()
