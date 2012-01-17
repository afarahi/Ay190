#!/usr/bin/python
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc
from NumIntegral import *
from math import *



n = 200;
xGridCount=zeros(n+1,int)
ErrSR=zeros(n+1)
ErrTR=zeros(n+1)

for i in range(0, n+1):
   xGridCount[i] = 10 + 10*i;

for j in range(0, n+1):
   x = linspace(0.0, pi, xGridCount[j]+1);
   fx = zeros(xGridCount[j]+1);
   for i in range(0, xGridCount[j]+1):
      fx[i] = x[i] * sin(x[i])
   IntegralSum = 0.0;
   for i in range(0, xGridCount[j]):
      IntegralSum = IntegralSum + TrapezoidalRule(fx[i],fx[i+1],x[i],x[i+1]);
   ErrTR[j] = abs(pi - IntegralSum) / pi

   IntegralSum = 0.0;
   i = 0
   while (i <= xGridCount[j]-1):
      IntegralSum = IntegralSum + SimpsonRule2(fx[i],fx[i+1],fx[i+2],x[i],x[i+2]);
      i += 2;
   ErrSR[j] = abs(pi - IntegralSum) / pi

figure(1)
semilogy(xGridCount,ErrSR,'r-',label='Simpson\'s Rule')
semilogy(xGridCount,ErrTR,'b-',label='Trapezoidal Rule')
xlabel('Grid Size')
ylabel('Relative Error')
legend(loc=1)
savefig('plot1.pdf')
