#!/opt/local/bin/python

import sys
from numpy import *
from math  import sin, pi, cos, exp
from scipy import *
from pylab import *
from matplotlib.pyplot import *

tmax = 4.0

def set_grid(tmax,nzones):
    t   = linspace(0.0,tmax,nzones+1);
    dt  = t[1] - t[0];
    return (t,dt)

def tov_BW(dt,olddata):
    newdata    = zeros(2)
    a          = (1.0-100.0*dt) 
    [Y1 , Y2]  = olddata
    newdata[0] = ((1.0-99.0*dt)*Y1 - 99.0*dt*Y2)/a
    newdata[1] = (-dt*Y1 + (1.0-dt)*Y2)/a
    return newdata

def tov_integrate(tmax,nzones):
    # set up grid
    (t,dt) = set_grid(tmax,nzones)
    
    # initialize some variables    
    tovdata = zeros((nzones,2))
    # 0 -- Y1
    # 1 -- Y2
 
    tovout  = zeros((nzones,3))
    # 0 -- Y1
    # 1 -- Y2
    # 2 -- t
    
    # initial conditions
    tovdata[0,0] = 1.0
    tovdata[0,1] = 0.0

    tovout[0,0] = tovdata[0,0]
    tovout[0,1] = tovdata[0,1]

    for i in range(nzones):
        tovout[i,2] = t[i]

    for i in range(nzones-1):
        tovdata[i+1,:] = tov_BW(dt,tovdata[i,:])
        
    tovout[:,0]  = tovdata[:,0]
    tovout[:,1]  = tovdata[:,1]

    return (tovout,dt)

# for convergence: 
# number of points
na = array([10000, 20000, 30000, 40000, 50000, 60000])#, 70000, 80000, 100000, 120000, 140000, 150000])
# to store masses
masses = zeros(len(na))
# to store the dts
dts = zeros(len(na))

####
### Exact Solution
####
n       = 10
(t,dt)  = set_grid(tmax,n)
Y1      = zeros(n+1)
Y2      = zeros(n+1) 
for i in range(0,n+1):
    Y1[i] = abs(exp(100.0*t[i]) / 100.0 + 99.0 / 100.0)
    Y2[i] = abs(- exp(100.0*t[i]) / 100.0 + 1.0 / 100.0)

figure(1)
semilogy(t,Y1,'o',label='Analytical Solution')
xlabel('t')
ylabel('Absolute Value of Y1')
figure(2)
semilogy(t,Y2,'o',label='Analytical Solution')
xlabel('t')
ylabel('Absolute Value of Y2')

RE = zeros(len(na))

for i in range(len(na)):
    print i+1
    n        = na[i]
    tov      = zeros((n,3))
    (tov,dt) = tov_integrate(4.0,n)
    figure(1)
    semilogy(tov[:,2] ,abs(tov[:,0]) ,label='Grid Size : %i' %(n))
    figure(2)
    semilogy(tov[:,2] ,abs(tov[:,1]) ,label='Grid Size : %i' %(n))
    RE[i]    = abs((tov[na[i]-1,0]-Y1[10])/Y1[10])


figure(1)
legend(loc=2)
savefig('plot4.pdf')
figure(2)
legend(loc=2)
savefig('plot5.pdf')
figure(3)
semilogy(na,RE)
xlabel('Gridd Size')
ylabel('Relative Error for explicit Euler integration')
savefig('plot6.pdf')

