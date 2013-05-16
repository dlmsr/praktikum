# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from linregress import linear_fit

# Laden der Geräte-Daten
L, C, R, R_ap = np.loadtxt('daten_und_Rap.txt', unpack=True)

########################################################################
# Teil a)

t, U = np.loadtxt('einhuellende.txt', unpack=True)
t *= 1e-6 # micro seconds
x = np.linspace(t.min(), t.max())
res = np.zeros((2, 4))

plt.plot(t, U, 'x')

for i in range(2):
    (m, b), err = linear_fit(t[i::2], np.log(U[i::2]))
    plt.plot(x, np.exp(m*x + b))

    res[i,0] = m
    res[i,1] = err[0]
    res[i,2] = b
    res[i,3] = err[1]

plt.yscale('log')
plt.grid(True, 'minor')
plt.grid(True, 'major')
plt.title('Bestimmung des Exponenten')
plt.ylabel('$U_C/\\mathrm{V}$')
plt.xlabel('$t/\\mathrm{s}$')
plt.savefig('exp-plot.pdf')
plt.close()

# Mittelwert und Standardabweichung vom Mittelwert der Steigungen
m_mean   = np.mean(res[:,0])
m_stderr = stats.sem(res[:,0])

# effektiver Dämpfungswiderstand
R_eff = 2*m_mean*L

print('Ergebnis der Ausgleichsrechnung:')
print('m, delta m, b, delta b')
print(res)
print('\nMittelwert: \t\t\t\t{0:0.4e}'.format(m_mean))
print('Standard-Abweichung: \t\t\t{0:0.4e}'.format(m_stderr))
print('effektiver Dämpfungswiderstand: \t{0:0.4e}'.format(R_eff))
print('Abweichung vom berechneten Wert: \t{0:0.4}'.format(abs(R - R_eff)/R))
print(72*'-')

########################################################################
# Teil b)


# Berechne R_ap aus Gerätedaten
R_ap_theo = 2*np.sqrt(L/C)

print('Berechneter R_ap: \t{0:0.4e}'.format(R_ap_theo))
print('Relative Abweichung R_ap gemessen/berechnet: {0:0.4f}'\
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

sec = np.where(U >= 65)
plt.subplot(212)
plt.plot(f[sec], U[sec]/U_0, '+')
plt.ylabel('$U/U_0$')
plt.xlabel('$\\nu/\\mathrm{Hz}$')
plt.grid(True, which='both')

plt.savefig('resonanz.pdf')
plt.close()

q_theo = 1/R*np.sqrt(L/C)
q = U.max()/U_0
breite_theo = R/L
breite = 2*np.pi*(f[sec].max() - f[sec].min())

print('Berechnete Resonanzüberhöhung: {0:0.4e}'.format(q_theo))
print('Experimentell bestimmte Resonanzüberhöhung: {0:0.4e}'.format(q))
print('Abweichung: {0:0.4f}'.format(abs(q_theo-q)/q_theo))
print('Berechnete Breite der Resonanzkurve: {0:0.4e}'.format(breite_theo))
print('Experimentell bestimmte Breite: {0:0.4e}'.format(breite))
print('Abweichung: {0:0.4f}'.format(abs(breite_theo-breite)/breite_theo))
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
