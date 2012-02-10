import sys
from pylab import *
from scipy import *
from numpy import *
from math  import exp, sqrt, pi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x0min = -1.0
x0max =  1.0
y0min = -1.0
y0max =  1.0
z0min = -1.0
z0max =  1.0

def load_data1():
   # Loading Data
   A   = loadtxt("write1.dat")
   x = list(A[:,0])
   y = list(A[:,1])
   z = list(A[:,2])
   return (x , y , z)

def load_data2():
   # Loading Data
   A   = loadtxt("write2.dat")
   x = list(A[:,0])
   y = list(A[:,1])
   z = list(A[:,2])
   return (x , y , z)


(Fx1 , Fy1 , Fz1)    = load_data1()
(Fx2 , Fy2 , Fz2)    = load_data2()


n = zeros(len(Fx1),int)
F1=zeros(len(Fx1))
F2=zeros(len(Fx1))
Er=zeros(len(Fx1))
print Fx1[0]
print Fx1[2]
for i in range(len(Fx1)):
    n[i] = i
    F1[i] = sqrt(Fx1[i]**2 + Fy1[i]**2 + Fz1[i]**2)
    F2[i] = sqrt(Fx2[i]**2 + Fy2[i]**2 + Fz2[i]**2)
    Er[i] = abs((F1[i] - F2[i])/F1[i])

figure()
plot(n,F1,'r-',n,F2,'b-')
show()
figure()
plot(n,Er,'r-')
xlabel('n (particle label)')
ylabel('relative error')
title('Critical angel = 0.8')
show()
#savefig('plot3.pdf')
