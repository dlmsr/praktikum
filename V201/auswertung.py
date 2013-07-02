import numpy as np
import matplotlib.pyplot as plt


def temperatur_in_kelvin(U):
    return 25.157*U - 0.19*U**2


# Lade Meßdaten
thermo_proben = np.loadtxt('thermo_proben.txt')
thermo_dewar  = np.loadtxt('thermo_dewar.txt')
daten_proben  = np.loadtxt('daten_proben.txt')

# Berechne Temperaturen 
temp_proben = temperatur_in_kelvin(thermo_proben[1:,:]).T
temp_dewar  = temperatur_in_kelvin(thermo_dewar[1:,:]).T

# Berechne Masse des verwendeten Wassers
V_w, c_w, rho_w = thermo_proben[1,:].T
V_x, V_y, V_g   = thermo_dewar[1,:].T

m_w = rho_w * V_w
m_x = rho_w * V_x
m_y = rho_w * V_y

m_k, rho_k, M_k, alpha, kappa = daten_proben.T

T_w, T_k, T_m = temp_proben
T_x, T_y, T_ms = temp_dewar

# Berechne Wärmekapazität und spez. Wärmkap. bei konstanten Druck
C_g = (c_w * m_y * (T_y - T_ms) - c_w * m_x * (T_ms - T_x))/(T_ms - T_x)
c_k = (c_w * m_w + C_g.mean()) * (T_m - T_w)/(T_k - T_m)/np.repeat(m_k, [3,1,1])

# Berechne Wärmekapazitäten bei konstantem Volumen
C_k_V = c_k*np.repeat(M_k, [3, 1, 1]) - np.repeat(9*alpha**2 * kappa * M_k/rho_k, [3, 1, 1]) * (T_m + 273.16)

print(C_g)
print(c_k)
print(C_k_V)
