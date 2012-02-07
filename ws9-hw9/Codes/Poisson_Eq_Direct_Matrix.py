import sys
from pylab import *
from numpy import *
from math  import exp, sqrt, pi
from matplotlib.pyplot import *
import scipy
import scipy.interpolate

#Constants
G   = 6.67300e-11

def Load_Data():
    # Loading Data
    A   = loadtxt("presupernova.dat")
    rho = A[0:500,4]
    r   = A[0:500,2]
    return [r , rho]

def Interpolation_Data(n):
    #Interpolation
    [r,rho]= Load_Data() 
    rep    = scipy.interpolate.splrep(r,rho,s=0)
    rnew   = linspace(0,10.0**9,n+1)
    dr     = rnew[1] - rnew[0]
    for i in range(n):
        rnew[i] = rnew[i] + dr*0.5
    rhonew = scipy.interpolate.splev(rnew,rep,der=0)
    for i in range(n):
        if (rnew[i] < r[0]):
            rhonew[i] = (rho[1]-rho[0])*(rnew[i]-r[0])/(r[1]-r[0]) + rho[0]
        else:
            break
    return [rnew,rhonew,dr]

def mass_cal(rho, r, dr):
    mass     = 0.0
    for i in range(len(r)):
        mass = mass + 4.0 * pi * r[i] * r[i] * dr * rho[i]
    return mass

def Set_Matrix(rho, r, dr):
    n           = len(r)
    A           = zeros([n,n])
    B           = zeros(n)
    mass        = mass_cal(rho, r, dr)
    for i in range(n-1):
        B[i]    = 4.0 * pi * G * rho[i]
    #Boundary Condition
    A[0,0]      = -1.0/(r[0]*dr) - 1.0/(dr*dr)
    A[0,1]      =  1.0/(r[0]*dr) + 1.0/(dr*dr)
    B[n-1]      = -G*mass/r[n-1]
    A[n-1,n-1]  =  1.0
    for i in range(1,n-1):
        A[i,i]  = -2.0/(dr*dr)
        A[i,i+1]=  1.0/(r[i]*dr) + 1.0/(dr*dr)
        A[i,i-1]= -1.0/(r[i]*dr) + 1.0/(dr*dr)
    return [A , B]

n          = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000]

[r,rho,dr] = Interpolation_Data(n[0])
[A , B]    = Set_Matrix(rho, r, dr)
phi        = linalg.solve(A,B)

for i in range(1,len(n)): 
   [r,rho,dr] = Interpolation_Data(n[i])
   [A , B]    = Set_Matrix(rho, r, dr)
   phinew     = linalg.solve(A,B)
   semilogy(n[i],abs((max(phinew)-max(phi))/max(phinew)),'bo')
   phi        = phinew

xlabel('n')
ylabel('Relative Error')
savefig('plot5.pdf')


figure()
semilogx(r,abs(phi))
savefig('plot3.pdf')


#################################
#################################
#################################
#################################
#################################
