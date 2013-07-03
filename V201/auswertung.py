import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import C2K, R
from scipy import stats

def temperatur_in_kelvin(U):
    return C2K(25.157*U - 0.19*U**2)

def rel_abw(theorie, exp):
    return np.abs(theorie - exp)/theorie

# Lade Meßdaten
thermo_proben = np.loadtxt('thermo_proben.txt')
thermo_dewar  = np.loadtxt('thermo_dewar.txt')
daten_proben  = np.loadtxt('daten_proben.txt')

# Berechne Temperaturen in Kelvin
temp_proben = temperatur_in_kelvin(thermo_proben[1:,:]).T
temp_dewar  = temperatur_in_kelvin(thermo_dewar[1:,:]).T

# Berechne Masse des verwendeten Wassers
V_w, c_w, rho_w = thermo_proben[0,:].T
V_x, V_y, V_g   = thermo_dewar[0,:].T

m_x, m_y, m_w = rho_w * np.array([V_x, V_y, V_w])

m_k, rho_k, M_k, alpha, kappa = daten_proben.T

T_w, T_k, T_m  = temp_proben
T_x, T_y, T_ms = temp_dewar

# Berechne Wärmekapazität und spez. Wärmkap. bei konstanten Druck
c_g_m_g = (c_w * m_y * (T_y - T_ms) - c_w * m_x * (T_ms - T_x)) \
          / (T_ms - T_x)
#c_g_m_g = 1/(T_ms - T_x) * (m_x + m_y) * c_w * (0.5*(T_x - T_y) - T_ms)
c_k = (c_w * m_w + c_g_m_g.mean()) / np.repeat(m_k, [3,1,1]) \
      * (T_m - T_w)/(T_k - T_m)

# Berechne Wärmekapazitäten bei konstantem Volumen
V_0 = M_k/rho_k
C_p = c_k * np.repeat(M_k, [3, 1, 1])
C_V = C_p - np.repeat(9*alpha**2 * kappa * V_0, [3, 1, 1]) * T_m

print('c_g m_g')
print(c_g_m_g)
print('c_k, C_p, C_V')
print(np.column_stack((c_k, C_p, C_V)))

np.savetxt('temperaturen', np.hstack((thermo_proben[1:,:], temp_proben.T)), fmt='%.4e', delimiter=' & ', newline = ' \\\\\n ')

A = np.hstack((thermo_dewar[1:], temp_dewar.T))
np.savetxt('temp2', A, fmt='%.4e', delimiter=' & ', newline = ' \\\\\n')

A = C_V[0:2]
print('mittelwert')
print(A.mean())
print('standardabweichung')
print(stats.sem(A))

B = np.array([A.mean()])
B = np.concatenate((B, C_V[3:]))

print(rel_abw(3*R, B))

