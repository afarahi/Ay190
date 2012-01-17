#!/usr/bin/python

###Integral Methods

def TrapezoidalRule(y1,y2,x1,x2):
   h1 = x2 - x1
   integral = (x2-x1)*(y1+y2)/2.0
   return integral


def SimpsonRule(y1,y2,y3,x1,x2):
   h1 = x2 - x1
   h2 = x3 - x2
   a = (2.0*h2**2 + h1*h2 - h1**2) / (6.0*h2)
   b = (h1 + h2)**3 / (6.0*h1*h2)
   c = (-h2**2 + h1*h2 - 2.0*h1**2) / (6.0*h2)
   integral = a*y3 + b*y2 + c*y1
   return integral

def SimpsonRule2(y1,y2,y3,x1,x3):
   h = x3 - x1
   integral = h * (y1 + 4.0*y2 + y3) / 6.0
   return integral

def GaussianQuadrature(y1,y2,a,b):
   integral = (b-a)*(y1+y2)/2.0 
   return integral

def GaussianQuadratureX1(a,b):
   from math import sqrt
   X1 = (b-a)/(2.0*sqrt(3.0)) + (a+b)/2.0
   return X1

def GaussianQuadratureX2(a,b):
   from math import sqrt
   X2 = - (b-a)/(2.0*sqrt(3.0)) + (a+b)/2.0
   return X2


