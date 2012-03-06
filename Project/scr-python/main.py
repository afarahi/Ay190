import sys
import time
from pylab import *
from scipy import *
from numpy import *
from N_Body_oct_tree_time     import build_mass_tree, force_calculator_OCT 
from N_Body_particle_particle import force_calculator_D
from N_Body_PM                import force_calculator_MP
######################
###################### main programm
######################

xmin    =-1.0
xmax    = 1.0
ymin    =-1.0
ymax    = 1.0
zmin    =-1.0
zmax    = 1.0

n = [100, 250 , 500, 750 , 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000]
elapsed1 = zeros(len(n))
elapsed2 = zeros(len(n))

random.seed(7)
#(x , y , z)      = load_data()
x = zeros(n[2])
y = zeros(n[2])
z = zeros(n[2])
for i in range(len(x)):
   x[i] = 2.0*random.random() - 1.0
   y[i] = 2.0*random.random() - 1.0
   z[i] = 2.0*random.random() - 1.0


for j in range(len(n)):
   '''
   random.seed(7)
   #(x , y , z)      = load_data()
   x = zeros(n[j])
   y = zeros(n[j])
   z = zeros(n[j])
   for i in range(len(x)):
      x[i] = 2.0*random.random() - 1.0
      y[i] = 2.0*random.random() - 1.0
      z[i] = 2.0*random.random() - 1.0
   '''

#   start             = time.clock()
#   (mass_tree,x,y,z) = build_mass_tree(x0min,x0max,y0min,y0max,z0min,z0max,n[j])
#   print "Building mass tree is done"
#   F                 = zeros((len(x),3))
#   for i in range(len(x)):
#      (F[i,0] , F[i,1] , F[i,2]) = force_calculator_OCT(mass_tree,x,y,z,i)
#   elapsed1[j]       = (time.clock() - start)

   nn          = (j+1)*2     
   start       = time.clock()
   FPM         = force_calculator_MP(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,nn)
   elapsed1[j] = (time.clock() - start)

#   start          = time.clock()
#   F = zeros((len(x),3))
#   for i in range(len(x)):
#       (F[i,0] , F[i,1] , F[i,2]) = force_calculator_D(x,y,z,i)
#   elapsed2[j]    = (time.clock() - start)

   print j, nn, elapsed1[j]#, elapsed2[j]

'''
figure()
plot(n,elapsed1,'-b',label="Oct-tree method")
plot(n,elapsed2,'-r',label="Direct method")
xlabel("Number of Particle")
ylabel("Time")
legend(loc=2)
savefig('plottime.pdf')

(mass_tree,x,y,z) = build_mass_tree(x0min,x0max,y0min,y0max,z0min,z0max,500)
print "Building mass tree is done"
F1                = zeros((len(x),3))
F2                = zeros((len(x),3))

for i in range(len(x)):
   (F1[i,0] , F1[i,1] , F1[i,2]) = force_calculator(mass_tree,x,y,z,i)

for i in range(len(x)):
    (F2[i,0] , F2[i,1] , F2[i,2]) = force_calculator_D(x,y,z,i)

print " done "

F11=zeros(len(x))
F21=zeros(len(x))
Er=zeros(len(x))

n = zeros(len(x),int)
for i in range(len(x)):
    n[i] = i
    F11[i] = sqrt(F1[i,0]**2 + F1[i,1]**2 + F1[i,2]**2)
    F21[i] = sqrt(F2[i,0]**2 + F2[i,1]**2 + F2[i,2]**2)
    Er[i] = abs((F11[i] - F21[i])/F11[i])

figure()
plot(n,F11,'r-',n,F21,'b-')
show()
figure()
plot(n,Er,'r-')
xlabel('n (particle label)')
ylabel('relative error')
title('Critical angel = 0.8')
show()'''
