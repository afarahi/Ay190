import sys
from pylab     import *
from numpy     import *
from math      import pi, cos, sin
from cmath     import exp

def Fourier_Trans(x,y,z,field,L,n): # L = max of deltaX or range of X, or Y, or Z
    n -= 0
    n3 = n**3
    l  = 0
    kx = zeros(n3) # kx
    ky = zeros(n3) # ky
    kz = zeros(n3) # kz
    cn = zeros(n3,complex) # Fourier Factor
    if (mod(n,2) == 0):
       for i in range(n):
          for j in range(n):
             for k in range(n):
                kx[l]   = - pi*n/L + 2.0*pi*i/L
                ky[l]   = - pi*n/L + 2.0*pi*j/L
                kz[l]   = - pi*n/L + 2.0*pi*k/L
                for m in range(len(x)):
                   cn[l] += exp(- 1.0j*kx[l]*x[m] - 1.0j*ky[l]*y[m] - 1.0j*kz[l]*z[m])*field[m]
                cn[l] = cn[l] / n3
                l += 1
    else:
       for i in range(n):
          for j in range(n):
             for k in range(n):
                kx[l]   = - pi*(n+1)/L + 2.0*pi*i/L
                ky[l]   = - pi*(n+1)/L + 2.0*pi*j/L
                kz[l]   = - pi*(n+1)/L + 2.0*pi*k/L
                for m in range(len(x)):
                   cn[l] += field[m]*exp(- 1.0j*kx[l]*x[m] - 1.0j*ky[l]*y[m] - 1.0j*kz[l]*z[m])
                cn[l] = cn[l] / n3
                l += 1
    return (cn, kx, ky, kz)

def Inv_Fourier_Trans(cn,kx,ky,kz,x,y,z):
    n     = len(cn)
    nn    = len(x)
    l     = 0
    field = zeros([nn],complex) # scalar field
    for i in range(nn):
       for m in range(n):
          field[i] += cn[m]*exp(1.0j*kx[m]*x[i] + 1.0j*ky[m]*y[i] + 1.0j*kz[m]*z[i])
    return field

def Inv_Fourier_Trans_dx(cn,kx,ky,kz,x,y,z):
    n     = len(cn)
    nn    = len(x)
    l     = 0
    field = zeros([nn],complex) # scalar field
    for i in range(nn):
       for m in range(n):
          field[i] += 1.0j*kx[m]*cn[m]*exp(1.0j*kx[m]*x[i] + 1.0j*ky[m]*y[i] + 1.0j*kz[m]*z[i])
    return field


def Fourier_Trans_inf(x,y,z,field,L,n): # L = max of deltaX or range of X, or Y, or Z
    n -= 0
    n3 = n**3
    l  = 0
    kx = zeros(n3) # kx
    ky = zeros(n3) # ky
    kz = zeros(n3) # kz
    cn = zeros(n3,complex) # Fourier Factor
    for i in range(n):
       for j in range(n):
          for k in range(n):
             kx[l]   = 2.0*pi*i/n
             ky[l]   = 2.0*pi*j/n
             kz[l]   = 2.0*pi*k/n
             for m in range(len(x)):
                cn[l] += exp(- 1.0j*kx[l]*x[m] - 1.0j*ky[l]*y[m] - 1.0j*kz[l]*z[m])*field[m]
             cn[l] = cn[l] * (L**3) / n3
             l += 1
    return (cn, kx, ky, kz)

def Inv_Fourier_Trans_inf(cn,kx,ky,kz,x,y,z):
    dk3   = (kx[1] - kx[0])**3
    n     = len(cn)
    nn    = len(x)
    l     = 0
    field = zeros([nn],complex) # scalar field
    for i in range(nn):
       for m in range(n):
          field[i] += cn[m]*exp(1.0j*kx[m]*x[i] + 1.0j*ky[m]*y[i] + 1.0j*kz[m]*z[i])*dk3
    return field

def Inv_Fourier_Trans_dx_inf(cn,kx,ky,kz,x,y,z):
    dk3   = (kx[1] - kx[0])**3
    n     = len(cn)
    nn    = len(x)
    l     = 0
    field = zeros([nn],complex) # scalar field
    for i in range(nn):
       for m in range(n):
          field[i] += 1.0j*kx[m]*cn[m]*exp(1.0j*kx[m]*x[i] + 1.0j*ky[m]*y[i] + 1.0j*kz[m]*z[i])*dk3
    return field


