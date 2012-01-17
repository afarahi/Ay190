#!/usr/bin/python
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc
from NumIntegral import *
from math import *
import scipy 
import scipy.special.orthogonal as scipy_orth 

########################################
########################################
#############  Constants  ##############
########################################
########################################

c = 299792458
M = 10**6
beta = 1.0/(20.0*M)
hbar = 6.58211928*(10**(-16)) # eV

########################################
########################################
####### Problem 2 - first Part #########
########################################
########################################

xGridCount = 50;

constant = 8.0*pi*((1.0/beta)**3) / ((2*pi*c*hbar)**3)

[laguerre_roots,laguerre_weights] = scipy_orth.l_roots(xGridCount,0)

IntegralSum = 0.0;

for i in range(0, xGridCount):
   IntegralSum = IntegralSum + laguerre_weights[i]*(laguerre_roots[i]**2 * exp(laguerre_roots[i]) / (exp(laguerre_roots[i]) + 1.0));
else:
   print "The Integral is (gaussian-laguerre method) : ", IntegralSum

##################################################
# Checking the integral answer by another method #
##################################################

xGridCount = 10000;

IntegralSum = 0.0;

x = linspace(0,50, xGridCount+1)
for i in range(0, xGridCount-1):
   x1 = x[i];
   x2 = x[i+1];
   y1 = (x[i]**2) / (exp(x[i]) + 1.0);
   y2 = (x[i+1]**2) / (exp(x[i+1]) + 1.0);
   IntegralSum = IntegralSum + TrapezoidalRule(y1,y2,x1,x2);
else:
   print "The Integral is (Trapezoidal Rule method) : ", IntegralSum

print "n_ = ", constant*IntegralSum

########################################
########################################
####### Problem 2 - Second Part ########
########################################
########################################

## Formula ::
## cp**2 dcp / (exp(beta cp) + 1)
## (E**2) / (exp (beta E) + 1)
## d(cp) = dE
##

xGridCount = 10;
EGridCount = 150;
Emax = 450*M

E = linspace(0,Emax, EGridCount+1)
fE= zeros(EGridCount+1)
dE = E[1] - E[0]

[legendre_roots,legendre_weights] = scipy_orth.p_roots(xGridCount,0)

constant = 8.0*pi/((2*pi*c*hbar)**3)

IntegralSum = 0.0;

for i in range(1, EGridCount):
   IntegralPart = 0.0;
   for j in range(0, xGridCount):
      X = E[i] - dE/2.0 + dE*legendre_roots[j]/2.0
      IntegralPart = IntegralPart + legendre_weights[j]*(X**2 / (exp(beta*X) + 1.0));
   IntegralSum = IntegralSum + dE*IntegralPart/2.0;
   fE[i] = IntegralPart/(2.0)
else:
      print "n_ = ", constant*IntegralSum


figure(1)
plot(E/M,fE,'b-')
xlabel('Energy (MeV)')
ylabel('Electron density spectral distribution')
savefig('plot2.pdf')
