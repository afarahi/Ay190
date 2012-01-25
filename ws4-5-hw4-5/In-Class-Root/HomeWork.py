#!/opt/local/bin/python

import sys
from numpy import *
from math  import sin, pi, cos
from scipy import *
from pylab import *
from matplotlib.pyplot import *

#Creating the poly. function
poly  = [0.0 , 0.0 , 0.0 , -1.0 , 5.0 , 3.0]

#calculating the nth deriviative
def fnp(pol,n):
    y = zeros(len(pol))
    for i in range(len(pol)):
        y[i] = pol[i]
    if (n > 0):
        for k in range(n): 
            for i in range(len(pol)-1):
                y[i] = (i+1)*y[i+1]
            y[len(pol)-1] = 0.0
    return y

def matfnp(pol):
    mat = zeros([len(pol)-1,len(pol)])
    for i in range(len(pol)-1):
        mat[i,:] = fnp(pol,i)
    return mat

def fp(pol,X):
    Y = 0.0
    for i in range(1,len(pol)):
        Y += i*pol[i]*X**(i-1)
    return Y

def f(pol,X):
    Y = 0.0
    for i in range(len(pol)):
        Y += pol[i]*X**i
    return Y

def NewtonMethod(pol,x0):
    i = 0
    x1 = x0
    x0 = 10.0
    while ( (abs((x1 - x0)/(x1+10e-30)) > 10e-20) and (i < 1000) ):
        x0 = x1
        i = i + 1
        x1 = x0 - f(pol,x0)/fp(pol,x0)
    return x1

def findrootII(pol,Xext1,Xext2):
    Yext1 = f(pol,Xext1)
    Yext2 = f(pol,Xext2)
    if (Yext1*Yext2 < 0e-10):
        x = NewtonMethod(pol,(Xext1+Xext2)/2.0)
    else:
        x = 1001.0
    return x

def findrootIp(pol,Xext):
    Yext = f(pol,Xext)
    Yinf = f(pol,10.0e10)
    if (Yext*Yinf < 1.0e-10):
        x = NewtonMethod(pol,Xext+0.1)
    else:
        x = 1001.0
    return x

def findrootIn(pol,Xext):
    Yext = f(pol,Xext) 
    Yinf = f(pol,-10.0e10)
    if (Yext*Yinf < 1.0e-10):  
        x = NewtonMethod(pol,Xext-0.1)
    else:
        x = 1001.0
    return x

def findroot0(pol):
    x = NewtonMethod(pol,-0.1)
    return x


Mat = matfnp(poly)

roots    = zeros(2)
roots[0] = findrootIp( Mat[len(poly)-3,:] , findroot0(Mat[len(poly)-2,:]) )
roots[1] = findrootIn( Mat[len(poly)-3,:] , findroot0(Mat[len(poly)-2,:]) ) 
roots.sort()
while(roots[len(roots)-1] == 1001.0):
    roots2 = roots[0:len(roots)-1]
    roots  = roots2 

for i in range(len(poly)-3):
    roots2 = zeros(len(roots)+1)
    if (len(roots) == 1):
        roots2[len(roots)] = findrootIp(Mat[len(poly)-4-i,:],roots[len(roots)-1])
        roots2[0]          = findrootIn(Mat[len(poly)-4-i,:],roots[0]           )
    else:
        for j in range(len(roots)-1):
            roots2[j+1]    = findrootII(Mat[len(poly)-4-i,:],roots[j],roots[j+1])
        roots2[len(roots)] = findrootIp(Mat[len(poly)-4-i,:],roots[len(roots)-1])
        roots2[0]          = findrootIn(Mat[len(poly)-4-i,:],roots[0]           )
    roots = roots2
    roots.sort()
    while(roots[len(roots)-1] == 1001.0):
        roots2 = roots[0:len(roots)-1]
        roots  = roots2 


print "Roots are : ", roots

