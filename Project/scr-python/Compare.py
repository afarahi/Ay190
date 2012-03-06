import sys
from pylab import *
from scipy import *
from numpy import *
from math  import exp, sqrt, pi

def load_data1():
   # Loading Data
   A   = loadtxt("time_MP.dat")
   n   = A[:,1]
   t1  = A[:,2]
   t2  = A[:,3]
   return (n,t1,t2)

def load_data2():
   # Loading Data
   A   = loadtxt("time_MPII.dat")
   n   = A[:,1]
   t   = A[:,2]
   return (n,t)

#(n,t1,t2) = load_data1()
(n,t1) = load_data2()
t2     = zeros(len(t1))
for i in range(len(t1)):
   n[i]  = n[i]**3
   t2[i] = 0.00075*n[i] * log(n[i]) 
#   n[i]  = n[i]**3

#plot(n,t1,'.b',label='Mesh Particle method')
#plot(n,t2,'.r',label='Direct Method')
plot(n,t1,'.b',label='Mesh Particle')
plot(n,t2,'.r',label='N log N')
xlabel('Number of Grid Cells')
ylabel('time')
legend(loc=2)
savefig('timeII.pdf')
show()


