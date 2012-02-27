#!/usr/bin/env python
import sys,math
from pylab import *
import numpy.random as random

random.seed(1)

n = 1000
m = 0
pi_cal = zeros(n)
ERR    = zeros(n)
for i in range(n):
   x = random.rand()
   y = random.rand()
   if (y <= sqrt(1 - x**2) ):
     m += 1
   pi_cal[i] = 4.0*m/(i+1)
   ERR[i]    = abs((pi_cal[i]-pi)/pi)

print pi_cal[n-1]
semilogy(ERR,'r-')
ylabel('Relative Error')
xlabel('n')
savefig('cal_pi.pdf')
show()
