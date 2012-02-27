#!/usr/bin/env python
import sys,math
from pylab import *
import numpy.random as random

random.seed(1)
n = 10000
L = 0.01
d = 1.0
x = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.5):
     x[i] = x[i-1] + L
   else:
     x[i] = x[i-1] - L
plot(x,'r')
x = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.5):
     x[i] = x[i-1] + L
   else:
     x[i] = x[i-1] - L
plot(x,'g')
x = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.5):
     x[i] = x[i-1] + L
   else:
     x[i] = x[i-1] - L
plot(x,'b')
x = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.5):
     x[i] = x[i-1] + L
   else:
     x[i] = x[i-1] - L
plot(x,'m')
savefig('random_walk_1D.pdf')
show()
