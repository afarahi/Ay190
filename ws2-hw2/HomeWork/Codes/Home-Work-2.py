#!/usr/bin/python
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc
from scipy.interpolate import interp1d

rc('text', usetex=True)
rc('font', family='serif')

mydata_all = loadtxt("data.dat",comments = '#')
tn = mydata_all[:,0]
amn = mydata_all[:,1]
tGridSize = size(tn)

n = 200;
t = linspace(0.0, 1.0, n+1);
am = zeros(n+1);
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
#########Lagrange Interpolation############
for i in range(0, n + 1):
   for j in range(0, tGridSize):
      L = 1.0;
      for k in range(0, tGridSize):
         if (j != k):
            L = L * (t[i] - tn[k]) / (tn[j] - tn[k]);
         else:
            L = L;
      am[i] = am[i] + L * amn[j];
else:
   print "Loop 1 is finished"	

figure(1)
plot(tn, amn, 'bo', label='Data' )
plot(t, am, 'r-', label='Lagrange Interpolation' )
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
############Linear Interpolation###########
for i in range(0, n + 1):
   for j in range(1, tGridSize):
      if (t[i] <= tn[j]):
         break
   am[i] = (t[i]-tn[j-1])/(tn[j]-tn[j-1])*amn[j] + (t[i]-tn[j])/(tn[j-1]-tn[j])*amn[j-1];
else:
   print "Loop 2 is finished"

figure(1)
plot(t, am, 'm-', label='Linear Interpolation')

figure(2)
plot(tn, amn, 'bo', label='Data' )
plot(t, am, 'm-', label='Linear Interpolation')
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
########Quadratic Interpolation############
for i in range(0, n + 1):
   for j in range(2, tGridSize):
      if (t[i] <= tn[j]):
         break
   a1 = (t[i]-tn[j-1])/(tn[j-2]-tn[j-1])*(t[i]-tn[j])/(tn[j-2]-tn[j]);
   a2 = (t[i]-tn[j-2])/(tn[j-1]-tn[j-2])*(t[i]-tn[j])/(tn[j-1]-tn[j]);
   a3 = (t[i]-tn[j-2])/(tn[j]-tn[j-2])*(t[i]-tn[j-1])/(tn[j]-tn[j-1]);
   am[i] = a1*amn[j-2] + a2*amn[j-1] + a3*amn[j];
else:
   print "Loop 3 is finished"

figure(1)
plot(t, am, 'k-', label='Quadratic Interpolation' )

figure(2)
plot(t, am, 'k-', label='Quadratic Interpolation' )

###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
#########Cubic Hermite Interpolation#######
ampn = zeros(tGridSize);

for i in range(0, tGridSize-1):
   ampn[i] = ( amn[i+1] - amn[i] ) / ( t[i+1] - t[i] ) 
else:
   ampn[tGridSize-1] = ampn[tGridSize-2];
   print "Forward method is used for finding the first derivative"
   print "Loop 4 is finished"

for i in range(0, n + 1):
   for j in range(1, tGridSize):
      if (t[i] <= tn[j]):
         break
   z = (t[i]-tn[j-1])/(tn[j]-tn[j-1])
   h00 = (2.0*(z**3) - 3.0*(z**2) + 1.0)
   h10 = (z**3 - 2.0*(z**2) + z)
   h01 = (-2.0*(z**3) + 3.0*(z**2))
   h11 = (z**3 - z**2) 
   am[i] = h00*amn[j-1] + h10*(tn[j]-tn[j-1])*ampn[j-1] + h01*amn[j] + h11*(tn[j]-tn[j-1])*ampn[j];
else:
   print "Loop 5 is finished"
figure(1)
plot(t, am, 'g-', label='Cubic Hermite Interpolation' )

###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
#########Cubic Hermite Interpolation#######
am = interp1d(tn, amn, kind='cubic')
figure(1)
plot(t, am(t), 'y-', label='Cubic Spline Interpolation' )
figure(2)
plot(t, am(t), 'r-', label='Cubic Spline Interpolation' )

figure(1)
xlabel(r'\textbf{time}')
ylabel(r'\textbf{apparent magnitude}')
legend(loc=2)
savefig('plot5.pdf')


figure(2)
xlabel(r'\textbf{time}')
ylabel(r'\textbf{apparent magnitude}')
legend(loc=2)
savefig('plot6.pdf')

