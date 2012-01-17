#!/usr/bin/python
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

n = 100;
tGridSize = 8;
tn = zeros(tGridSize+2);
amn = zeros(tGridSize+2);
t = linspace(0.0, 1.0, n+1);
am = zeros(n+1);


tn[0] = 0.0;
tn[1] = 0.2;
tn[2] = 0.3;
tn[3] = 0.4;
tn[4] = 0.5;
tn[5] = 0.6;
tn[6] = 0.7;
tn[7] = 0.8;
tn[8] = 1.0;
tn[9] = 1.2;

amn[0]= 0.302;
amn[1]= 0.185;
amn[2]= 0.106;
amn[3]= 0.093;
amn[4]= 0.240;
amn[5]= 0.579;
amn[6]= 0.561;
amn[7]= 0.468;
amn[8]= 0.302;
amn[9]= 0.185;

#mydata_all = loadtxt("cop.dat",command = "#")
#time = mydata_all[:,0]
#time = mydata_all[:,1]

for i in range(0, n + 1):
   for j in range(0, tGridSize + 1):
      L = 1.0;
      for k in range(0, tGridSize + 1):
         if (j != k):
            L = L * (t[i] - tn[k]) / (tn[j] - tn[k]);
         else:
            L = L;
      am[i] = am[i] + L * amn[j];
else:
   print "Loop 1 is finished"	


figure(1)
plot(tn, amn, 'bo', label='Data' )
plot(t, am, 'r-', label='Polynomial Interpolation' )
xlabel(r'\textbf{time}')
ylabel(r'\textit{apparent magnitude}',fontsize=16)


for i in range(0, n + 1):
   for j in range(1, tGridSize + 1):
      if (t[i] <= tn[j]):
         break
   am[i] = (t[i]-tn[j-1])/(tn[j]-tn[j-1])*amn[j] + (t[i]-tn[j])/(tn[j-1]-tn[j])*amn[j-1];
else:
   print "Loop 2 is finished"

plot(t, am, 'm-', label='Linear Piecewise Interpolation' )

ampn = zeros(tGridSize+1);

for i in range(0, tGridSize):
   ampn[i] = ( amn[i+1] - amn[i] ) / ( t[i+1] - t[i] ) 
else:
   ampn[tGridSize] = ampn[tGridSize-1]
   print "Loop 3 is finished"


for i in range(0, n + 1):
   for j in range(1, tGridSize + 1):
      if (t[i] <= tn[j]):
         break
   z = (t[i]-tn[j-1])/(tn[j]-tn[j-1])
   s01 = (2.0*(z**3) - 3.0*(z**2) + 1.0)
   s02 = (2.0*((1.0-z)**3) - 3.0*((1.0-z)**2) + 1.0)
   s11 = (z**3 - 2.0*(z**2) + z)
   s12 = ((1.0-z)**3 - 2.0*((1.0-z)**2) + (1.0-z)) 
   am[i] = s01*amn[j-1] + s02*amn[j] + s11*(tn[j]-tn[j-1])*ampn[j-1] - s12*(tn[j]-tn[j-1])*ampn[j];
else:
   print "Loop 4 is finished"

plot(t, am, 'g-', label='Cubic Piecewise Interpolation' )



for i in range(0, n + 1):
   for j in range(1, tGridSize + 1):
      if (t[i] <= tn[j]):
         break
   a1 = (t[i]-tn[j])/(tn[j-1]-tn[j])*(t[i]-tn[j+1])/(tn[j-1]-tn[j+1]);
   a2 = (t[i]-tn[j-1])/(tn[j]-tn[j-1])*(t[i]-tn[j+1])/(tn[j]-tn[j+1]);
   a3 = (t[i]-tn[j-1])/(tn[j+1]-tn[j-1])*(t[i]-tn[j])/(tn[j+1]-tn[j]);
   am[i] = a1*amn[j-1] + a2*amn[j] + a3*amn[j+1];
else:
   print "Loop 4 is finished"

plot(t, am, 'k-', label='Quadratic Interpolation' )


legend(loc=2)
savefig('plot5.pdf')
