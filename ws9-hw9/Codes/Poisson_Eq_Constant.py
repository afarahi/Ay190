import sys
from pylab import *
from math  import exp, sqrt, pi
import scipy
import scipy.interpolate

#Constants
G   = 6.67300e-11

def rh_z_rep(z,rho,r,dr):
#    znew   = z + dr*(4.0*pi*G*rho - 2.0*z/r)
    znew   = (z + dr*4.0*pi*G*rho)/(1+dr*2.0/r)
    return znew

def rh_phi_rep(z,phi,dr):
    phinew = phi - dr*z
    return phinew

def mass_cal(rho, r, dr):
    mass = 0
    for i in range(len(r)):
        mass = mass + 4.0 * pi * r[i] * r[i] * dr * rho[i]
    return mass

def fw_z_sol(r,rho,dr):
    n      = len(r)
    z      = zeros(n)
    #initial condition
    z[0]   = 0.0
    for i in range(1,n):
        z[i] = rh_z_rep(z[i-1],rho[i],r[i],dr)
    return z

def fw_phi_sol(r,rho,dr):
    n        = len(r)
    phi      = zeros(n)
    #initial condition
    z        = fw_z_sol(r,rho,dr)
    mass     = mass_cal(rho, r, dr)
    print mass , len(z), len(r), len(phi)
    phi[n-1] = - G * mass / r[n-1] 
    for i in range(1,n):
        phi[n-1-i] = rh_phi_rep(z[n-i-1],phi[n-i],dr)
    return phi

n      = 100000
r      = linspace(0,10.0**9,n+1)
dr     = r[1] - r[0]
rho    = zeros(n+1)
for i in range(n+1):
    rho[i] = 5.0e8
#semilogx(r,rho)
#show()
phi = fw_phi_sol(r, rho, dr)

figure()
semilogy(r,abs(phi),'r-')
phi1 = zeros(n+1)
for i in range(n+1):
    phi1[i] = 2.0 * pi * G * rho[i] * (r[i]*r[i] - 3.0 * 10.0**9 * 10.0**9 ) / 3.0  
semilogy(r,abs(phi1),'b-')
show()

semilogy(r,abs((phi1-phi)/phi))
xlabel('R')
ylabel('Relative Error')
savefig('Plot1.pdf')
#################################
#################################
#################################
#################################
#################################
