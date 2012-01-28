#!/opt/local/bin/python

import sys
from numpy import *
from math  import sin, pi, cos
from scipy import *
from pylab import *
from matplotlib.pyplot import *
import time

data  = loadtxt("myburn_0002_z1.dat", comments='*')

figure(1)
semilogy(data[:,1],data[:,2])
xlabel('Time (s)')
ylabel('H1 mass fraction')
savefig('plot1.pdf')


figure(2)
semilogx(data[:,1],data[:,5])
xlabel('Time (s)')
ylabel('He4 mass fraction')
savefig('plot2.pdf')
