#!/usr/bin/python
from matplotlib.pyplot import *
from numpy import *
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

xGridSize = 100;
x = zeros(xGridSize+1);
fx = zeros(xGridSize+1);
dxF = zeros(xGridSize+1);
dxC = zeros(xGridSize+1);
dxR = zeros(xGridSize+1);
h = (6.0+2.0)/xGridSize;

for i in range(0, xGridSize + 1):
   x[i] = -2.0 + i*h;
   fx[i]= x[i]**3 - 5.0*(x[i]**2) + x[i];
else:
   print "Loop 1 is finished"


for i in range(1, xGridSize):
   dxF[i] = (fx[i+1] - fx[i]) / h;
   dxC[i] = (fx[i+1] - fx[i-1]) / (2.0*h);
   dxR[i] = 3.0*x[i]**2 - 10.0*x[i] + 1.0;
else:
   print "Loop 2 is finished"

figure(1)
plot(x[1:xGridSize-1], abs(dxC[1:xGridSize-1]-dxR[1:xGridSize-1]), 'r', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{Absolute error}',fontsize=16)
title(r"Absolute error of first derivative (central differencing method)", fontsize=16, color='b')
grid(True)


figure(2)
plot(x[1:xGridSize-1], abs(dxF[1:xGridSize-1]-dxR[1:xGridSize-1]), 'r', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{Absolute error}',fontsize=16)
title(r"Absolute error of first derivative (forward differencing method)", fontsize=16)
grid(True)

figure(3)
plot(x[1:xGridSize-1], dxC[1:xGridSize-1], 'r', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (central differencing method)", fontsize=16)
grid(True)


figure(4)
plot(x[1:xGridSize-1], dxF[1:xGridSize-1], 'r', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (forward differencing method)", fontsize=16)
grid(True)


xGridSize = 200;
x = zeros(xGridSize+1);
fx = zeros(xGridSize+1);
dxF = zeros(xGridSize+1);
dxC = zeros(xGridSize+1);
dxR = zeros(xGridSize+1);
h = (6.0+2.0)/xGridSize;

for i in range(0, xGridSize + 1):
   x[i] = -2.0 + i*h;
   fx[i]= x[i]**3 - 5.0*(x[i]**2) + x[i];
else:
   print "Loop 1 is finished"


for i in range(1, xGridSize):
   dxF[i] = (fx[i+1] - fx[i]) / h;
   dxC[i] = (fx[i+1] - fx[i-1]) / (2.0*h);
   dxR[i] = 3.0*x[i]**2 - 10.0*x[i] + 1.0;
else:
   print "Loop 2 is finished"


figure(1)
plot(x[1:xGridSize-1], abs(dxC[1:xGridSize-1]-dxR[1:xGridSize-1]), 'b', label='Grid Size = %i' %(xGridSize))
xlabel(r'\textbf{x}')
ylabel(r'\textit{Absolute error}',fontsize=16)
grid(True)
legend()
savefig('plot1.pdf')

figure(2)
plot(x[1:xGridSize-1], abs(dxF[1:xGridSize-1]-dxR[1:xGridSize-1]), 'b', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{Absolute error}',fontsize=16)
grid(True)
legend()
savefig('plot2.pdf')


figure(3)
plot(x[1:xGridSize-1], dxC[1:xGridSize-1], 'b', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (central differencing method)", fontsize=16, color='b')
grid(True)


figure(4)
plot(x[1:xGridSize-1], dxF[1:xGridSize-1], 'b', label='Grid Size = %i' %(xGridSize) )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (central differencing method)", fontsize=16)
grid(True)


figure(3)
plot(x[1:xGridSize-1], dxR[1:xGridSize-1], 'g', label='Analytical' )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (Analytical)", fontsize=16, color='k')
grid(True)
legend()
savefig('plot3.pdf')

figure(4)
plot(x[1:xGridSize-1], dxR[1:xGridSize-1], 'g', label='Analytical' )
xlabel(r'\textbf{x}')
ylabel(r'\textit{df/dy}',fontsize=16)
title(r"First derivative (Analytical)", fontsize=16)
legend()
savefig('plot4.pdf')
