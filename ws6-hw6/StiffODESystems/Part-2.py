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

def tov_RHS(t,data):
    rhs       = zeros(2)
    [Y1 , Y2] = data
    rhs[0]    = Y1 - 99.0 * Y2
    rhs[1]    = - Y1 + 99.0 * Y2
    return rhs

def tov_RK4(old_data,t,dt):
    k1          = zeros(2)
    k2          = zeros(2)
    k3          = zeros(2)
    k4          = zeros(2)
    new_data    = zeros(2)
    [Y1 , Y2]   = tov_RHS(t,old_data)
    k1[0]       = dt*Y1;
    k1[1]       = dt*Y2;
    [Y1 , Y2]   = tov_RHS(t+0.5*dt,old_data+0.5*k1)
    k2[0]       = dt*Y1;
    k2[1]       = dt*Y2;
    [Y1 , Y2]   = tov_RHS(t+0.5*dt,old_data+0.5*k2)
    k3[0]       = dt*Y1;
    k3[1]       = dt*Y2;
    [Y1 , Y2]   = tov_RHS(t+dt,old_data+k3)
    k4[0]       = dt*Y1;
    k4[1]       = dt*Y2;
    new_data[0] = old_data[0] + (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
    new_data[1] = old_data[1] + (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;
    return new_data


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
        tovdata[i+1,:] = tov_RK4(tovdata[i,:],t[i],dt)
        tovout[i+1,0] = tovdata[i+1,0]
        tovout[i+1,1] = tovdata[i+1,1]

    return (tovout,dt)

# for convergence: 
# number of points
na = array([5, 10, 50, 100, 200, 500, 1000, 5000])#, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000])
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
    n = na[i]
    tov = zeros((n,3))
    (tov,dt) = tov_integrate(4.0,n)
    figure(1)
    semilogy(tov[:,2] ,abs(tov[:,0]) ,label='Grid Size : %i' %(n))
    figure(2)
    semilogy(tov[:,2] ,abs(tov[:,1]) ,label='Grid Size : %i' %(n))
    RE[i] = abs((tov[na[i]-1,0]-Y1[10])/Y1[10])


figure(1)
legend(loc=2)
savefig('plot1.pdf')
figure(2)
legend(loc=2)
savefig('plot2.pdf')
figure(3)
semilogy(na,RE)
xlabel('Gridd Size')
ylabel('Relative Error for Backward Euler method')
savefig('plot3.pdf')

