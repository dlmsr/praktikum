# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt

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

I_S = np.array(I_S)
print(I_S)

########################################################################
# Teil b)

U, I = UI[-1].T
plt.figure()
plt.plot(U, I**(2/3), 'x')

plt.close()
