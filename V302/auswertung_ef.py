from numpy import *
from matplotlib.pyplot import *
import linregress

v,U_br = loadtxt("wien_robinson.txt",unpack = True)

U_s = 2 #Volt
U_br = U_br *10**(-3) #Volt
v = 2*pi*v #Kreisgeschwindigkeit

U = U_br/U_s

v0 = 1/(1000*294.75*10**(-9))

omega = v/v0

semilogx(omega, U, '*', label = "Messwerte")

def theo(om):
    return sqrt((om**2 -1)**2 /(9*( (1-om**2)**2 + 9*om**2)))

x = linspace(0,55,10000)
semilogx(x, theo(x), label = "Theorie")
#Bem.: Messwerte erreichen nie 0, da Oberwellen auftreten!

xlabel("$\omega$/$\omega_0$")
ylabel("$U_{Br}/U_S$")
legend(loc = "lower right")
grid()
savefig("wien_robinson_plot.pdf")
close()

#Klirrfaktor k

U0= 1.4 *10**(-3)
k = 1/(U_s/(sqrt(45)*U0) -1)

