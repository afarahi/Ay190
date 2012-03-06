import sys
from pylab import *
from scipy import *
from numpy import *
from math  import sqrt

def load_data():
   # Loading Data
   A = loadtxt("Data1.dat")
   x = list(A[:,0])
   y = list(A[:,1])
   z = list(A[:,2])
   #Processing Data
   i = 0
   while (i < len(x)):
      if ( (x[i]>1.0) or (x[i]<-1.0)  or (y[i]>1.0) or (y[i]<-1.0) or (z[i]>1.0) or (z[i]<-1.0) ):
         x.pop(i)
         y.pop(i)
         z.pop(i)
      else:
         i+=1
   return (x , y , z)

def force_calculator_D(x,y,z,i):
    G      = 6.67e-11 # Gravitational constant
    mass   = 1        # Mass of each particle
    #Initializing
    Fx     = 0.0
    Fy     = 0.0
    Fz     = 0.0
    for j in range(0,i):
       r   = sqrt((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2)
       Fx += G*mass*mass*(x[j] - x[i])/(r**3)
       Fy += G*mass*mass*(y[j] - y[i])/(r**3)
       Fz += G*mass*mass*(z[j] - z[i])/(r**3)
    for j in range(i+1,len(x)):
       r   = sqrt((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2)
       Fx += G*mass*mass*(x[j] - x[i])/(r**3)
       Fy += G*mass*mass*(y[j] - y[i])/(r**3)
       Fz += G*mass*mass*(z[j] - z[i])/(r**3)
    return (Fx,Fy,Fz)
