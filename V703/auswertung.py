# encoding: utf-8
import numpy as np
from linregress import linear_fit
import matplotlib.pyplot as plt

# Lade Meßdaten
messung_teil_a_d = np.loadtxt('teil_a_und_d.txt')
messung_teil_b_c = np.loadtxt('teil_b_und_c.txt')

teil_b, zweiquellen, osci = np.vsplit(messung_teil_b_c, np.array([1,4]))
Z, U, I = messung_teil_a_d.T


fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)

plt.xlim(200, 800)
plt.suptitle('Charakteristik des Zählrohrs', fontsize=16)

ax1.grid()

ax1.plot(U, Z, '.-')

# zweite Hälfte der Meßwerte herunterziehen
Z = np.concatenate((Z[:8], Z[8:] - 300))
ax1.plot(U[7:], Z[7:], 'r.-')

ax2.plot(U[:2], Z[:2], 'xb')
ax2.plot(U[2:12], Z[2:12], 'xr')
ax2.plot(U[12:], Z[12:], 'xb')

x = U[2:12]
y = Z[2:12]
m, s_m, b, s_b = linear_fit(x, y)

x = np.linspace(300, 700)
ax2.plot(x, m*x+b, '--g')
ax2.grid()
plt.xlim(200, 800)

print('Länge des Plateau-Bereichs: {:.3f}'.format(U[12]-U[3]))
print('Plateau-Steigung in % pro 100 V: {:.3%}'.format(m*100))

# Teil b)

print('Zeitl. Abstand Primär-Nachentladungen: {:.3e}'.format(teil_b[0,0]*teil_b[0,1]))

# Teil c)
# Berechne Totzeiten

DT = zweiquellen[0,1]
N1, N2, N12 = zweiquellen[:,0]/DT

T = (N1 + N2 - N12)/(2*N1*N2) # Totzeit in Sekunden

print('Zwei Quellen-Methode')
print('Totzeit in Sekunden: {:.3e}'.format(T))
print('Oscilloscope')
print('Totzeit in Sekounden: {:.3e}'.format(osci[0,0]*osci[0,1]))
