import numpy as np
import matplotlib.pyplot as plt

messung_teil_a_d = np.loadtxt('teil_a_und_d.txt')

Z, U, I = messung_teil_a_d.T

plt.plot(U, Z, '+b')
plt.plot(U, Z, 'b')

plt.plot(U[8:], Z[8:] - 300, 'gx')
plt.plot(U[8:], Z[8:] - 300, 'g')

plt.grid()
plt.show()

