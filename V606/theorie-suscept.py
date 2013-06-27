# This work is licensed under the Creative Commons
# Attribution-NonCommercial 3.0 Unported License.  To view a copy of
# this license, visit http://creativecommons.org/licenses/by-nc/3.0/.

import numpy as np

from scipy.constants import mu_0, N_A, hbar, k, physical_constants

mu_B = physical_constants['Bohr magneton'][0]

proben = np.loadtxt('verbindungen.txt')

# Raumtemperatur in Kelvin
T = 297.16

res_matrix = []

for (J, L, S, rho, M) in proben:

    # Berechne Landé-Faktor
    g_J = 1/2 * (3 + (S*(S + 1) - L*(L + 1))/J/(J + 1))

    # Berechne Teilchendichte in 1/m^3
    N   = 2 * rho/M * N_A * 1e6

    # Berechne Suszeptibilität
    chi = mu_0 * mu_B**2 * g_J**2 * N * J * (J + 1)/(3 * k * T)

    res_matrix.append([rho, M, N, g_J, chi])

res_matrix = np.array(res_matrix)
print(res_matrix)
np.savetxt('theorie-suscept.txt', res_matrix, delimiter=" & ",
           newline="\\\\\n", fmt="%.4e")
