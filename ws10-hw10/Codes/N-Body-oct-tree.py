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

def load_data():
   # Loading Data
   A   = loadtxt("Data1.dat")
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

def mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z):
   n = 0
   xc= 0.0
   yc= 0.0
   zc= 0.0
   for i in range(len(x)):
      if ((x[i]>=xmin) and (x[i]<xmax)):
         if((y[i]>=ymin) and (y[i]<ymax)):
            if((z[i]>=zmin) and (z[i]<zmax)):
               n += 1
               xc+= x[i]
               yc+= y[i]
               zc+= z[i]
   if (n != 0):
      xc = xc/n
      yc = yc/n
      zc = zc/n
   else:
      xc = 0.0
      yc = 0.0
      zc = 0.0
   return (n, xc, yc, zc)

def ini_tree(x1min,x1max,y1min,y1max,z1min,z1max,x,y,z):
    (n,xc,yc,zc)= mass_center_and_counter(x1min,x1max,y1min,y1max,z1min,z1max,x,y,z)
    mass_tree   = {'data':(x1min,x1max,y1min,y1max,z1min,z1max,n,xc,yc,zc)}
    if (n > 1):
       #Creat_Childrens
       #1st Children
       xmin        = x1min
       ymin        = y1min
       zmin        = z1min
       xmax        = x1min + (x1max - x1min)/2.0
       ymax        = y1min + (y1max - y1min)/2.0
       zmax        = z1min + (z1max - z1min)/2.0
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({1:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #2nd Children
       xmin        = x1min + (x1max - x1min)/2.0
       xmax        = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({2:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #3rd Children
       xmin        = x1min
       xmax        = x1min + (x1max - x1min)/2.0
       ymin        = y1min + (y1max - y1min)/2.0
       ymax        = y1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({3:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #4th Children
       xmin        = x1min + (x1max - x1min)/2.0
       xmax        = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({4:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #5th Children
       xmin        = x1min
       ymin        = y1min
       zmin        = z1min + (z1max - z1min)/2.0
       xmax        = x1min + (x1max - x1min)/2.0
       ymax        = y1min + (y1max - y1min)/2.0
       zmax        = z1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({5:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #6th Children
       xmin        = x1min + (x1max - x1min)/2.0
       xmax        = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({6:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #7th Children
       xmin        = x1min
       xmax        = x1min + (x1max - x1min)/2.0
       ymin        = y1min + (y1max - y1min)/2.0
       ymax        = y1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({7:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
       #8th Children
       xmin        = x1min + (x1max - x1min)/2.0
       xmax        = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       mass_tree.update({8:{'data':(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc)}})
    return mass_tree

def creat_mass_children(mass_tree,path,x1min,x1max,y1min,y1max,z1min,z1max,x,y,z):
       #Creat_Childrens
       #1st Children
       xmin        = x1min
       ymin        = y1min
       zmin        = z1min
       xmax        = x1min + (x1max - x1min)/2.0
       ymax        = y1min + (y1max - y1min)/2.0
       zmax        = z1min + (z1max - z1min)/2.0
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({1:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #2nd Children
       xmin = x1min + (x1max - x1min)/2.0
       xmax = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({2:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #3rd Children
       xmin = x1min
       xmax = x1min + (x1max - x1min)/2.0
       ymin = y1min + (y1max - y1min)/2.0
       ymax = y1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({3:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #4th Children
       xmin = x1min + (x1max - x1min)/2.0
       xmax = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({4:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #5th Children
       xmin = x1min
       ymin = y1min
       zmin = z1min + (z1max - z1min)/2.0
       xmax = x1min + (x1max - x1min)/2.0
       ymax = y1min + (y1max - y1min)/2.0
       zmax = z1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({5:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #6th Children
       xmin = x1min + (x1max - x1min)/2.0
       xmax = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({6:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #7th Children
       xmin = x1min
       xmax = x1min + (x1max - x1min)/2.0
       ymin = y1min + (y1max - y1min)/2.0
       ymax = y1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({7:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       #8th Children
       xmin = x1min + (x1max - x1min)/2.0
       xmax = x1max
       (n,xc,yc,zc)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
       s    = "mass_tree" + path + ".update({8:{'data':(" + str(xmin) + "," + str(xmax) + "," + str(ymin) + "," + str(ymax) + "," + str(zmin) + "," + str(zmax) + "," + str(n) + "," + str(xc) + "," + str(yc) + "," + str(zc) + ")}})"
       exec(s)
       return mass_tree

def add_index(path):
    path     += "[1]"
    return path

def continue_index(path):
    n            = len(path)
    if (int(path[0-2]) < 8):
       path      = path[0:n-3] + "[" + str(int(path[n-2])+1) + "]"
    elif(path == "[8]"):
       path      = "END"  
    else:
       path      = continue_index(path[0:n-3])
    return path

def build_mass_tree(xmin,xmax,ymin,ymax,zmin,zmax):
    #Loading Data
    (x , y , z) = load_data()
    print "data reading done"
    mass_tree   = ini_tree(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
    #Initializing Data
    mass_num    = 0
    print "Initializing done"
    path    = "[1]"
    b       = 1
    while(b == 1):
       s    = "n = mass_tree" + path + "['data'][6]"
       exec(s)
       if (n > 1):
          s = "x1min = mass_tree" + path + "['data'][0]"
          exec(s)
          s = "x1max = mass_tree" + path + "['data'][1]"
          exec(s)
          s = "y1min = mass_tree" + path + "['data'][2]"
          exec(s)
          s = "y1max = mass_tree" + path + "['data'][3]"
          exec(s)
          s = "z1min = mass_tree" + path + "['data'][4]"
          exec(s)
          s = "z1max = mass_tree" + path + "['data'][5]"
          exec(s)
          mass_tree = creat_mass_children(mass_tree,path,x1min,x1max,y1min,y1max,z1min,z1max,x,y,z)
          path      = add_index(path)
       elif (n == 1):
          path      = continue_index(path)
          mass_num += 1
          print mass_num
       else:
          path      = continue_index(path)
       if (path == "END"):
           b = 0
    return (mass_tree, mass_num)

def find_path(mass_tree,x,y,z):
    path = "[1]"
    b = 1
    while (b == 1):
       s          = "n = mass_tree" + path + "['data'][6]"
       exec(s)
       if (n > 1):
          path    = add_index(path)
       elif (n == 1):
          s       = "xc = mass_tree" + path + "['data'][7]"
          exec(s)
          s       = "yc = mass_tree" + path + "['data'][8]"
          exec(s)
          s       = "zc = mass_tree" + path + "['data'][9]"
          exec(s)
          if (x == xc and y == yc and z == zc):
             b    = 0
          else:
             path = continue_index(path)
       else:
          path    = continue_index(path)
       if (path == "END"):
           b      = 0
    return path

def force_calculator(mass_tree,x,y,z,i):
    G      = 6.67e-11 # Gravitational constant
    mass   = 1.0      # Mass of each particle
    angel  = 0.8      # Critical Angel
    #Finding the path of particle i
    path_i = find_path(mass_tree,x[i],y[i],z[i])
    #Initializing
    Fx     = 0.0
    Fy     = 0.0
    Fz     = 0.0
    path   = "[1]"
    b      = 1
    while(b == 1):
       s          = "xc = mass_tree" + path + "['data'][7]"
       exec(s)
       s          = "yc = mass_tree" + path + "['data'][8]"
       exec(s)
       s          = "zc = mass_tree" + path + "['data'][9]"
       exec(s)
       s          = "mass_num = mass_tree" + path + "['data'][6]"
       exec(s)
       r          = sqrt((xc - x[i])**2 + (yc - y[i])**2 + (zc - z[i])**2)
       if ((r > 0.0) and (mass_num != 0)):
          s       = "d = mass_tree" + path + "['data'][1]" + " - mass_tree" + path + "['data'][0]"
          exec(s)
          theta   = d/r
          if ((theta < angel) or (mass_num == 1)):
             Fx  += G*mass*mass*mass_num*(xc - x[i])/(r**3)
             Fy  += G*mass*mass*mass_num*(yc - y[i])/(r**3)
             Fz  += G*mass*mass*mass_num*(zc - z[i])/(r**3)
             path = continue_index(path)
          else:
             path = add_index(path)   
       else:
          path    = continue_index(path)
       if (path == "END"):
           b = 0
    return (Fx,Fy,Fz)


######################
###################### main programm
######################

(x , y , z)    = load_data()
(mass_tree, n) = build_mass_tree(x0min,x0max,y0min,y0max,z0min,z0max)

F = zeros((len(x),3))
for i in range(len(x)):
    (F[i,0] , F[i,1] , F[i,2]) = force_calculator(mass_tree,x,y,z,i)
#    print i , " done"

savetxt('write2.dat',F)
