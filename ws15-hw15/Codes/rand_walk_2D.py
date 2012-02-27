#!/usr/bin/env python
import sys,math
from pylab import *
import numpy.random as random

#random.seed(1)
n = 10000
L = 0.01
d = 1.0
x = zeros(n)
y = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.125):
     x[i] = x[i-1] + L
     y[i] = y[i-1]
   elif (ran < 0.25):
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.375):
     x[i] = x[i-1]
     y[i] = y[i-1] + L
   elif (ran < 0.5):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.625):
     x[i] = x[i-1] - L
     y[i] = y[i-1]
   elif (ran < 0.75):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
   elif (ran < 0.875):
     x[i] = x[i-1]
     y[i] = y[i-1] - L
   else:
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
plot(x,y,'b')
x = zeros(n)
y = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.125):
     x[i] = x[i-1] + L
     y[i] = y[i-1]
   elif (ran < 0.25):
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.375):
     x[i] = x[i-1]
     y[i] = y[i-1] + L
   elif (ran < 0.5):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.625):
     x[i] = x[i-1] - L
     y[i] = y[i-1]
   elif (ran < 0.75):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
   elif (ran < 0.875):
     x[i] = x[i-1]
     y[i] = y[i-1] - L
   else:
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
plot(x,y,'r')
x = zeros(n)
y = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.125):
     x[i] = x[i-1] + L
     y[i] = y[i-1]
   elif (ran < 0.25):
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.375):
     x[i] = x[i-1]
     y[i] = y[i-1] + L
   elif (ran < 0.5):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.625):
     x[i] = x[i-1] - L
     y[i] = y[i-1]
   elif (ran < 0.75):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
   elif (ran < 0.875):
     x[i] = x[i-1]
     y[i] = y[i-1] - L
   else:
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
plot(x,y,'g')

x = zeros(n)
y = zeros(n)
for i in range(1,n):
   ran = random.rand()
   if (ran < 0.125):
     x[i] = x[i-1] + L
     y[i] = y[i-1]
   elif (ran < 0.25):
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.375):
     x[i] = x[i-1]
     y[i] = y[i-1] + L
   elif (ran < 0.5):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] + L/sqrt(2)
   elif (ran < 0.625):
     x[i] = x[i-1] - L
     y[i] = y[i-1]
   elif (ran < 0.75):
     x[i] = x[i-1] - L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
   elif (ran < 0.875):
     x[i] = x[i-1]
     y[i] = y[i-1] - L
   else:
     x[i] = x[i-1] + L/sqrt(2)
     y[i] = y[i-1] - L/sqrt(2)
plot(x,y,'m')
savefig('random_walk_2D.pdf')
show()
