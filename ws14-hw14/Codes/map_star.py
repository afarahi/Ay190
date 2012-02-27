import sys,math
from pylab    import *

rhomin = 1.0e6
ggrav  = 6.672e-8
dt     = 1.0e-5
K1     = 1.2435e15 * (0.5e0**(4.0/3.0))
gamma1 = 1.28


def setup_star(hyd):
    poly = loadtxt('poly.dat')
    prad = poly[:,0]
    prho = poly[:,1]
    nn = len(prho)
    for i in range(hyd.n):
        hyd.rho[i] = max(rhomin,linterp(hyd.x[i],nn,prho,prad))
        hyd.press[i] = K1 * hyd.rho[i]**gamma1
        hyd.eps[i] = hyd.press[i] / (gamma1-1.0) / hyd.rho[i]
        if(hyd.rho[i] <= rhomin):
            hyd.rho[i] = hyd.rho[i] / 5.0    
    return hyd


def linterp(xx,n,f,x):
    i=0
    while(i<n and x[i] < xx):
        i=i+1
    if(i==n):
        ff = rhomin
    elif(i==0):
        ff = (f[1]-f[0])/(x[1]-x[0]) * (xx - x[0]) + f[0]
    else:
        ff = (f[i]-f[i-1])/(x[i]-x[i-1]) * (xx-x[i-1]) + f[i-1]
    return ff


