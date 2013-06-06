# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scs
from scipy.constants import k, e, sigma, m_e, pi, hbar

########################################################################
# Teil a)

# Kennlinien laden
kennlinien = np.loadtxt('kennlinien.txt')

UI = np.hsplit(kennlinien, 5)
I_S = []

for (U, I) in map(np.transpose, UI):
    # nur die Werte ungleich 0
    U = U[U.nonzero()]
    I = I[I.nonzero()]
    # Sättigungsstrom speichern
    I_S.append(I[1:].max())

    plt.title('Heizspannung: {:.1f} Volt, '
              'Heizstrom: {:.1f} Ampère'.format(U[0], I[0]))
    plt.xlabel('Spannung in Volt')
    plt.ylabel('Strom in Ampère')
    plt.grid()
    plt.ylim(0, I[1:].max()*1.05)
    plt.plot(U[1:], I[1:], 'x')
    plt.axhline(I[1:].max(), linestyle='--')
    plt.savefig('kennlinie-' + '{:.1f}'.format(I[0]).replace('.', '_') +
                'A.pdf')
    plt.close()

I_S  = np.array(I_S)
U, I = np.vstack(np.split(kennlinien[0,:], 5)).T
tab  = np.array([U, I, I_S]).T
np.savetxt('saettig.txt', tab, delimiter=' & ', newline=' \\\\\n',
           fmt='%.4e')

########################################################################
# Teil b)

# Meßwerte der max. Heizleistung
# ohne Heizspannung/strom
U, I = UI[-1].T
U = U[1:]
I = I[1:]

# lineare Regression
x = np.log(U[:12])
y = np.log(I[:12])

m, b, r, prob2, see = scs.linregress(x, y)

# Standardabweichung der Parameter
mx   = x.mean()
sx2  = ((x - mx)**2).sum()
sd_m = see * np.sqrt(1./len(x) + mx*mx/sx2)
sd_b = see * np.sqrt(1./sx2)

print('Raumladungsgesetz:')
print('m = {:.3f} +- {:.3f}'.format(m, sd_m))
print('b = {:.3f} +- {:.3f}'.format(b, sd_b))
print('|m - 3/2|/(3/2) = {:.2%}'.format(abs(m-1.5)/1.5))

plt.loglog(U, I, '+')
plt.plot(np.exp(x), np.exp(m*x+b))
plt.grid(which='both')
plt.title('Auffindung des Raumladungsgebiets')
plt.xlabel('Spannung in Volt')
plt.ylabel('Strom in Ampère')
plt.savefig('raumladung.pdf')
plt.close()

########################################################################
# Teil c)

UI   = np.loadtxt('anlaufmessung.txt')
U, I = UI.T

# lineare Regression
x = -U
y = np.log(I)

m, b, r, prob2, see = scs.linregress(x, y)

# Standardabweichung der Parameter
mx   = x.mean()
sx2  = ((x - mx)**2).sum()
sd_m = see * np.sqrt(1./len(x) + mx*mx/sx2)
sd_b = see * np.sqrt(1./sx2)

plt.semilogy(-U, I, '+')
plt.semilogy(-U, np.exp(m*x+b))
plt.grid(which='both')
plt.savefig('anlaufgebiet.pdf')
plt.close()

DeltaT = e/(m**2*k)*sd_m
print('Kathodentemperatur: {:.3f} +- {:.3f}K'.format(e/(m*k), DeltaT))

########################################################################
# Teil d)
eta = 0.28
f = 0.35e-4
UI = np.loadtxt('kennlinien.txt')

# Heizspannungen und -ströme
U, I = np.vstack(np.split(UI[0,:], 5)).T

T = ((U*I - 1)/(f*eta*sigma))**(1/4)

np.savetxt('temp.txt', np.array([U, I, T]).T, delimiter=' & ',
           newline=' \\\\\n', fmt='%.4e')

#########################################################################
# Teil e)
A = 4*pi*e*m_e*k**2/hbar**3
austritt = -k*T*np.log(2*I_S/(f*A*T**2))/e

tab = np.array([I_S, T, austritt]).T
tab = np.append(tab, [tab.mean(axis=0), scs.sem(tab, axis=0)], axis=0)

np.savetxt('austritt.txt', tab, delimiter=' & ', newline=' \\\\\n',
           fmt='%.4e')
