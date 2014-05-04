# -*- encoding: utf-8 -*-

from __future__ import unicode_literals

from numpy import *
from matplotlib.pyplot import *

#Messwerte#

mikro = array([2.14,2.54,2.86,3.26,3.58,3.96])
swr=array([0,2,4,6,8,10])
swr = swr+9
eich=array([9,11.5,15,20,23.5,29])

plot(mikro,swr,'b*',label='SWR-Anzeige')
plot(mikro,eich,'r',label='Eichkurve')

xlabel('Mikrometerstellung / mm')
ylabel('Dämpfung / dB')
legend(loc='upper left')
grid()
savefig('daempfung.pdf')
