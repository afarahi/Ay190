#!/opt/local/bin/python

import sys
from numpy import *
from math  import sin, pi, cos
from scipy import *
from pylab import *
from matplotlib.pyplot import *
from GE import gauss
import time

start = time.clock()

mymatrixA5 = loadtxt("dataA5.dat")
mymatrixB5 = loadtxt("dataB5.dat")

elapsed = (time.clock() - start)

print elapsed

start = time.clock()

X = gauss(mymatrixA5,mymatrixB5)

elapsed = (time.clock() - start)

print elapsed
print X

start = time.clock()

X = linalg.solve(mymatrixA5,mymatrixB5)

elapsed = (time.clock() - start)

print elapsed
print X
