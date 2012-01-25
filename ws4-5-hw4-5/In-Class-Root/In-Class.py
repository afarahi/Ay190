#!/opt/local/bin/python

import sys
from numpy import *
from math  import sin, pi, cos
from scipy import *
from pylab import *
from matplotlib.pyplot import *

# global constants
period = 365.25635
omega  = 2.0 * pi / period
k      = 1000.0
a      = 1.496e6*k
e      = 0.0167
tp     = 1.0
t      = [91.0*tp , 182.0*tp , 273.0*tp]

def f(E,time,ep):
    y = E - omega * time - ep*sin(E)    
    return y

def fp(E,time,ep):
    yp = 1.0 - ep*cos(E)    
    return yp


for j in range(len(t)):
    i = 0
    E0 = 10.0
    E1 = 1.0
    while ( (abs((E1 - E0)/(E1+10e-30)) > 10e-10) and (i < 1000) ):
        E0 = E1
        i = i + 1
        E1 = E0 - f(E0,t[j],e)/fp(E0,t[j],e)
        figure(j)
        semilogy(i,abs((E1 - E0)/(E1+10e-30)),'ob')
        xlabel('n')
        ylabel('Relative Error')
    print 'time = ', t[j], ' E = ', E1, ' i = ', i
    savefig('plotI%i.pdf' %(j))

#part b 

print 'Part b'
print '---------------------'

e = 0.99999

for j in range(len(t)):
    i = 0
    E0 = 10.0
    E1 = 1.0
    while ( (abs((E1 - E0)/(E1+10e-30)) > 10e-10) and (i < 1000) ):
        E0 = E1
        i = i + 1
        E1 = E0 - f(E0,t[j],e)/fp(E0,t[j],e)
        figure(j+3)
        semilogy(i,abs((E1 - E0)/(E1+10e-30)),'ob')
        xlabel('n')
        ylabel('Relative Error')
    print 'time = ', t[j], ' E = ', E1, ' i = ', i
    savefig('plotII%i.pdf' %(j))

