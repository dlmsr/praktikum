# encoding: utf-8
import numpy as np
from linregress import linear_fit
from scipy.constants import e
import matplotlib.pyplot as plt

# Lade Meßdaten
messung_teil_a_d = np.loadtxt('teil_a_und_d.txt')
messung_teil_b_c = np.loadtxt('teil_b_und_c.txt')

teil_b, zweiquellen, osci = np.vsplit(messung_teil_b_c, np.array([1,4]))
Z, U, I = messung_teil_a_d.T


fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
plt.suptitle('Charakteristik des Zählrohrs', fontsize=16)

ax1.grid()
ax1.set_xlim(200, 800)
ax1.set_xlabel('Spannung in Volt')
ax1.set_ylabel('Anzahl der Impulse')

ax2.grid()
ax2.set_xlim(200, 800)
ax2.set_xlabel('Spannung in Volt')

# tatsächliche Meßwerte plotten
ax1.plot(U, Z, '.-')

# zweite Hälfte der Meßwerte herunterziehen
Z = np.concatenate((Z[:8], Z[8:] - 300))
ax1.plot(U[7:], Z[7:], 'r.-')

# Meßwerte für die lineare Ausgleichsrechnung markieren
# Fehlerbalken drangeben

# Spannung, Zählrate und stat. Fehler der Zählrate in eine Matrix
plot_data = np.column_stack((U, Z, np.sqrt(Z)))
split_list = np.hsplit(plot_data.T, [2, 12])

for x, y, yerr in split_list:
    if len(x) > 2:
        ax2.errorbar(x, y, yerr=yerr, fmt='xr')
    else:
        ax2.errorbar(x, y, yerr=yerr, fmt='xb')

# lineare Ausgleichsrechnung
x = U[2:12]
y = Z[2:12]
M = linear_fit(x, y)

m = M[0]
b = M[2]

x = np.linspace(250, 750)
ax2.plot(x, m*x+b, '--g')

plt.savefig('charakteristik.pdf')
plt.close()

print('Länge des Plateau-Bereichs: {:.3f}'.format(U[12]-U[3]))
print('Parameter der Ausgleichsgeraden')
print(M)
print('Plateau-Anstieg: {:.3f}'.format( m*100/Z[2]))

np.savetxt('teil_a.txt', np.array([U, Z, np.sqrt(Z)]).T, fmt='%d',
           delimiter=' & ', newline='\\\\\n')

# Teil b)

print('Zeitl. Abstand Primär-Nachentladungen: {:.3e}'.format(teil_b[0,0]
                                                             * teil_b[0,1]))

# Teil c)
# Berechne Totzeiten

DT = zweiquellen[0,1]
N1, N2, N12 = zweiquellen[:,0]/DT

T = (N1 + N2 - N12)/(2*N1*N2) # Totzeit in Sekunden

print('Zwei Quellen-Methode')
print('Totzeit in Sekunden: {:.3e}'.format(T))
print('Oscilloscope')
print('Totzeit in Sekunden: {:.3e}'.format(osci[0,0]*osci[0,1]))

# Teil d)

# Strom in microAmp und Ablesefehler des Stromes
I  = I * 1e-6
DI = 0.02e-6

# Zählrate
N  = Z/120
DN = np.sqrt(Z)/120

# Ladung pro Teilchen
Q = I/N
DQ = np.sqrt(( DI / N )**2 + (I/N**2 * DN)**2)

plt.plot(U, Q/e, '+')
plt.title('Pro Teilchen freigesetzte Ladung')
plt.xlabel('Zählrohrspannung U in Volt')
plt.ylabel('Freigesetzte Ladung in Elementarladungen')
plt.grid()
plt.savefig('freigesetzte_ladung.pdf')

np.savetxt('teil_d.txt', np.array([U, I*1e9, Q/e*1e-6, DQ/e*1e-6]).T, fmt='%d', 
           delimiter=' & ', newline='\\\\\n')
plt.close()
