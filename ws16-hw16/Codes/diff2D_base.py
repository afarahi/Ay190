#!/usr/bin/env python
#from __future__ import division
from matplotlib.patches import Patch
from pylab import *

def gauss(x,y,x0,y0,sig):
    Z = zeros((len(x),len(y)))
    for i in range(len(x)):
       for j in range(len(y)):
           Z[i,j] = exp(-((x[i]-x0)**2 + (y[j]-y0)**2) / (100*sig))
#    print "implement me"
    return Z

nx = 100
ny = 100
x0 = 50
y0 = 50
sig = 10.0
nit = 1000
D = 1.0

x = linspace(0, 100, nx)
y = linspace(0, 100, ny)
X,Y = meshgrid(x, y)

dx = x[1]-x[0]
dx2 = dx*dx
idx2 = 1.0/dx2
dt = 0.5 * dx**2 / (2.0 * D) #"implement me"

Z = gauss(x, y,x0,y0,sig)
RHS = zeros((nx,ny))

def calc_rhs(Z,nx,ny,dx2):
    RHS = zeros((nx,ny))
#    print "Implement me"
    for i in range(1,nx-1):
       for j in range(1,ny-1):
           RHS[i,j] =  D*(Z[i+1,j]-2.0*Z[i,j]+Z[i-1,j])/dx2 + D*(Z[i,j+1]-2.0*Z[i,j]+Z[i,j-1])/dx2 
    return RHS

def set_bound(Z,nx,ny):
    # boundaries
    Z[0,0] = Z[1,1]
    Z[0,:] = Z[1,:]
    Z[:,0] = Z[:,1]
    Z[nx-1,ny-1] = Z[nx-2,nx-1]
    Z[nx-1,:] = Z[nx-2,:]
    Z[:,ny-1] = Z[:,ny-2]
    return Z

for it in range(nit):
    RHS    = calc_rhs(Z,nx,ny,dx2)

    Z[:,:] = Z[:,:] + D*RHS*dt
    Z      = set_bound(Z,nx,ny)

    print it
    if it % 25 == 0:
        print "plotting"
        clf()
        pcolor(X, Y, Z, vmin=0.0, vmax=1.0)
        colorbar()
        fname = 'frame%04d.png'%(it/25)
        print 'Saving frame', fname
        savefig(fname)





