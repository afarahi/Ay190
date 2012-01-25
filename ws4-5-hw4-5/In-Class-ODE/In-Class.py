#!/opt/local/bin/python

import sys
from numpy import *
from scipy import *
from pylab import *
from matplotlib.pyplot import *

# global constants
ggrav = 6.67e-8
clite = 3.0e10
msun = 1.99e33

# EOS
# neutron stars:
# polyG = 2.0
# polyK = 100.0 * 5.55e38/6.1755e17**polyG

# EOS for 
# white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

# central values
rhoc = 1.0e10

# minimum pressure
rhomin = 10e-10 * rhoc
min_press = polyK * rhomin**polyG 

# grid
rmax = 1.0e9

def set_grid(rmax,nzones):
    # set up the grid and return the
    # radius array and dr
    rad = linspace(0.0,rmax,nzones+1);
    dr  = rad[1] - rad[0];
    return (rad,dr)

def tov_RHS(r,data):
    rhs    = zeros(2)
    mass   = data[1]
    press  = max(data[0],min_press)
    rho    = (press/polyK)**(1.0/polyG)
    if (abs(r) > 1.0):
        rhs[0] = - ggrav * mass * rho / (r**2)
    else:
        rhs[0] = 0.0
    rhs[1] = 4.0 * pi * rho * (r**2)
    return rhs

def tov_RK2(old_data,r,dr):
    k1          = zeros(2)
    k2          = zeros(2)
    new_data    = zeros(2)
    (p , m)     = tov_RHS(r,old_data)
    k1[0]       = dr*p;
    k1[1]       = dr*m;
    (p , m)     = tov_RHS(r+0.5*dr,old_data+0.5*k1)
    k2[0]       = dr*p;
    k2[1]       = dr*m;
    new_data[0] = old_data[0] + k2[0];
    new_data[1] = old_data[1] + k2[1];
    return new_data
    
def tov_RK3(old_data,r,dr):
    k1          = zeros(2)
    k2          = zeros(2)
    k3          = zeros(2)
    new_data    = zeros(2)
    [p , m]     = tov_RHS(r,old_data)
    k1[0]       = dr*p;
    k1[1]       = dr*m;
    [p , m]     = tov_RHS(r+0.5*dr,old_data+0.5*k1)
    k2[0]       = dr*p;
    k2[1]       = dr*m;
    [p , m]     = tov_RHS(r+dr,old_data-k1+2.0*k2)
    k3[0]       = dr*p;
    k3[1]       = dr*m;
    new_data[0] = old_data[0] + (k1[0] + 4.0*k2[0] + k3[0])/6.0;
    new_data[1] = old_data[1] + (k1[1] + 4.0*k2[1] + k3[1])/6.0;
    return new_data

def tov_RK4(old_data,r,dr):
    k1          = zeros(2)
    k2          = zeros(2)
    k3          = zeros(2)
    k4          = zeros(2)
    new_data    = zeros(2)
    [p , m]     = tov_RHS(r,old_data)
    k1[0]       = dr*p;
    k1[1]       = dr*m;
    [p , m]     = tov_RHS(r+0.5*dr,old_data+0.5*k1)
    k2[0]       = dr*p;
    k2[1]       = dr*m;
    [p , m]     = tov_RHS(r+0.5*dr,old_data+0.5*k2)
    k3[0]       = dr*p;
    k3[1]       = dr*m;
    [p , m]     = tov_RHS(r+dr,old_data+k3)
    k4[0]       = dr*p;
    k4[1]       = dr*m;
    new_data[0] = old_data[0] + (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
    new_data[1] = old_data[1] + (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;
    return new_data

def tov_integrate(rmax,nzones,k):
    # set up grid
    (rad,dr) = set_grid(rmax,nzones)
    
    # initialize some variables    
    tovdata = zeros((nzones,2))
    # 0 -- press
    # 1 -- mbary
 
    tovout  = zeros((nzones,5))
    # 0 -- rho
    # 1 -- press
    # 2 -- eps
    # 3 -- mass
    # 4 -- rad
     
    # central values
    tovdata[0,0] = polyK * rhoc**polyG
    tovdata[0,1] = 0.0

    for i in range(nzones-1):
        tovout[i,4] = rad[i]

    # you will need to track the surface (where press <= press_min)
    isurf = 0
    for i in range(nzones-1):
        # integrate one step using RK2 (can change this to
        # RK3 or RK4)
        if(k == 2):
            tovdata[i+1,:] = tov_RK2(tovdata[i,:],rad[i],dr)
        if(k == 3):
            tovdata[i+1,:] = tov_RK3(tovdata[i,:],rad[i],dr)
        if(k == 4):
            tovdata[i+1,:] = tov_RK4(tovdata[i,:],rad[i],dr)

        # check if press below 0
        if((tovdata[i+1,0] <= min_press) and (isurf == 0)):
            isurf = i

        # press and mass
        tovout[i+1,1] = tovdata[i+1,0]

        if (i+1 > isurf and isurf > 0):
           tovout[i+1,3] = tovdata[isurf,1]
        else:
           tovout[i+1,3] = tovdata[i+1,1]
        # compute density
        tovout[i+1,0] = (tovdata[i+1,0]/polyK)**(1.0/polyG)
        # compute eps
        tovout[i+1,2] = tovdata[i+1,1]/((polyG-1.0)*tovout[i+1,0])

    return (tovout,isurf,dr)

# for convergence: 
# number of points
na = array([  5, 10, 20, 50, 100, 500, 1000, 1500, 2000, 2500 , 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000])
# to store masses
masses = zeros(len(na))
# to store the drs
drs = zeros(len(na))

isurf_c = zeros(len(na))

for j in range(3):
    for i in range(len(na)):
        (tov_star,isurf,dr) = tov_integrate(rmax,na[i],j+2)
        isurf_c[i] = isurf 
#        print 'mass : ', tov_star[isurf,3]/msun, 'radios : ' , tov_star[isurf,4]/1000, ' Km', 'radios : ' , isurf*dr/1000.0, ' Km'
    for i in range(2, len(na)): 
        figure(3)
        if(j == 0):
           plot(na[i],abs((isurf_c[i]-isurf_c[i-1])/(isurf_c[i-1]-isurf_c[i-2])),'ro',label='K2')
        if(j == 1):
           plot(na[i],abs((isurf_c[i]-isurf_c[i-1])/(isurf_c[i-1]-isurf_c[i-2])),'bo',label='K3')
        if(j == 2):
           plot(na[i],abs((isurf_c[i]-isurf_c[i-1])/(isurf_c[i-1]-isurf_c[i-2])),'mo',label='K4')

    for i in range(1, len(na)):
        figure(4)
        if(j == 0):
           semilogy(na[i],abs((tov_star[isurf_c[i],3]-tov_star[isurf_c[i-1],3])/tov_star[isurf_c[i],3]),'ro',label='K2')
        if(j == 1):
           semilogy(na[i],abs((tov_star[isurf_c[i],3]-tov_star[isurf_c[i-1],3])/tov_star[isurf_c[i],3]),'bo',label='K3')
        if(j == 2):
           semilogy(na[i],abs((tov_star[isurf_c[i],3]-tov_star[isurf_c[i-1],3])/tov_star[isurf_c[i],3]),'mo',label='K4')


figure(3)
xlabel('grid size')
ylabel(' self-convergence factor')
savefig('plot3.pdf')
figure(4)
xlabel('grid size')
ylabel('relative error (For mass of the star)')
savefig('plot4.pdf')
figure(1)
subplot(2, 1, 1)
plot(tov_star[0:isurf,4] ,tov_star[0:isurf,3]/msun )
xlabel('r')
ylabel('mass / mass(sun)')
subplot(2, 1, 2)
semilogy(tov_star[0:isurf,4] ,tov_star[0:isurf,0])
xlabel('r')
ylabel('rho')
savefig('plot1.pdf')
figure(2)
semilogy(tov_star[0:isurf,4] ,tov_star[0:isurf,1])
xlabel('r')
ylabel('Pressure')
savefig('plot2.pdf')
#figure(3)
#semilogy(tov_star[:,4] ,tov_star[:,2])
#xlabel('r')
#ylabel('Energy')
