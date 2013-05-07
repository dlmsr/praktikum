# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so
from linregress import linear_fit

# Temp, Schwarz, Weiss, Metall matt, Metall glänzend
data = np.loadtxt("stefan-boltzmann.txt", unpack=True)

T0 = 293.16 # kelvin
T = data[0] + 273.16 # kelvin
U_therm = data[1:] * 1e-6 # volts

m = np.empty(4)
Dm = np.empty(4)
b = np.empty(4)
Db = np.empty(4)

epsilon = np.empty(4)

for i in range(4):

    (a, c), (a_err, c_err) = linear_fit(T**4-T0**4, U_therm[i])

    m[i] = a
    Dm[i] = a_err
    b[i] = c
    Db[i] = c_err
    epsilon[i] = a/m[0]

# Gausssche Fehlerformel
Depsilon = np.sqrt((Dm/m[0])**2 + (m*Dm[0]/m[0]**2)**2)

name = ["Schwarz", "Weiß", "Metall (matt)", "Metall (glänzend)"]

for i in range(4):
    print(name[i]+"\n")
    print("\tm = \t{0:0.4e} +- {1:0.4e}".format(m[i], Dm[i]))
    print("\tb = \t{0:0.4e} +- {1:0.4e}".format(b[i], Db[i]))
    print("\tepsilon = {0:0.4e} +- {1:0.4e}\n".format(epsilon[i], Depsilon[i]))

x = T**4 - T0**4

for i in range(4):
    plt.subplot(221+i)
    plt.xlabel("$(T^4 - T_0^4)/\mathrm{K}^4$")
    plt.ylabel("$U/\mathrm{V}$")
    plt.ticklabel_format(style='sci', scilimits=(-3,4), axis='both')
    plt.grid()

    plt.plot(x, U_therm[i], "+")
    plt.plot(x, m[i]*x + b[i])

plt.tight_layout()
plt.savefig("stefan-boltzmann.pdf")
plt.close()

data = np.loadtxt("abstandsmessung.txt", unpack=True)

x = (data[0] + 20) * 1e-2 # m
U = data[1:] * 1e-6 # volts


def f(x, a, b):
    return a/x**2 + b

for i in range(2):
    (a, b), cov = so.curve_fit(f, x, U[i])

    print("a = {0:0.4e}".format(a))
    print("b = {0:0.4e}".format(b))
    print("Da = {0:0.4e}".format(np.sqrt(cov[0][0])))
    print("Da/a = {0:0.4f}".format(np.sqrt(cov[0][0])/a))
    print("Db = {0:0.4e}".format(np.sqrt(cov[1][1])))

    plt.subplot(211+i)
    plt.plot(x, U[i], '+')

    y = np.linspace(0.1, 0.8)
    plt.plot(y, f(y, a, b))

    plt.xlim(0.1, 0.8)
    plt.grid()
    plt.ticklabel_format(style='sci', scilimits=(-3,4), axis='both')
    plt.xlabel("$x/\mathrm{m}$")
    plt.ylabel("$U/\mathrm{V}$")

plt.tight_layout()
plt.savefig("abstandsmessung.pdf")
plt.close()
