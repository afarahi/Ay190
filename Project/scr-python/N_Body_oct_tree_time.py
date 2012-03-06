import sys
from pylab import *
from scipy import *
from numpy import *
from math  import exp, sqrt, pi

def put_in_order_data(x,y,z):
   for i in range(0,len(x)):
      for j in range(0,len(x)-1):
         if (x[j] > x[j+1]):
            xp = x[j]
            yp = y[j]
            zp = z[j]
            x[j] = x[j+1]
            y[j] = y[j+1]
            z[j] = z[j+1]      
            x[j+1] = xp
            y[j+1] = yp
            z[j+1] = zp      
   return (x , y , z)

def load_data(n):
   # Loading Data
   A   = loadtxt("Data1.dat")
   x = list(A[0:n,0])
   y = list(A[0:n,1])
   z = list(A[0:n,2])
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

def mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index0=0):
   n     = 0
   xc    = 0.0
   yc    = 0.0
   zc    = 0.0
   index = index0
   for i in range(index0,len(x)):
      if (x[i]>=xmin):
         if (x[i]>=xmax):
            break;
         else:
            if((y[i]>=ymin) and (y[i]<ymax)):
               if((z[i]>=zmin) and (z[i]<zmax)):
                  n += 1
                  xc+= x[i]
                  yc+= y[i]
                  zc+= z[i]
      else:
         index += 1
   if (n != 0):
      xc = xc/n
      yc = yc/n
      zc = zc/n
   else:
      xc = 0.0
      yc = 0.0
      zc = 0.0
   return (n, xc, yc, zc , index)

def ini_tree(x1min,x1max,y1min,y1max,z1min,z1max,x,y,z):
     mass_tree      = {}
     #Creat_Childrens
     #1st Children
     xmin        = x1min
     ymin        = y1min
     zmin        = z1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymax        = y1min + (y1max - y1min)/2.0
     zmax        = z1min + (z1max - z1min)/2.0
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({1:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #2nd Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({2:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #3rd Children
     xmin        = x1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymin        = y1min + (y1max - y1min)/2.0
     ymax        = y1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({3:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #4th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({4:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #5th Children
     xmin        = x1min
     ymin        = y1min
     zmin        = z1min + (z1max - z1min)/2.0
     xmax        = x1min + (x1max - x1min)/2.0
     ymax        = y1min + (y1max - y1min)/2.0
     zmax        = z1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({5:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #6th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({6:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #7th Children
     xmin        = x1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymin        = y1min + (y1max - y1min)/2.0
     ymax        = y1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({7:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #8th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
     mass_tree.update({8:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     return mass_tree

def creat_mass_children(mass_tree,path,x,y,z):
     x1min= mass_tree[path][0]
     x1max= mass_tree[path][1]
     y1min= mass_tree[path][2]
     y1max= mass_tree[path][3]
     z1min= mass_tree[path][4]
     z1max= mass_tree[path][5]
     index= mass_tree[path][10]
     path = path*10
     #Creat_Childrens
     #1st Children
     xmin        = x1min
     ymin        = y1min
     zmin        = z1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymax        = y1min + (y1max - y1min)/2.0
     zmax        = z1min + (z1max - z1min)/2.0
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+1:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #2nd Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+2:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #3rd Children
     xmin        = x1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymin        = y1min + (y1max - y1min)/2.0
     ymax        = y1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+3:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #4th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+4:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #5th Children
     xmin        = x1min
     ymin        = y1min
     zmin        = z1min + (z1max - z1min)/2.0
     xmax        = x1min + (x1max - x1min)/2.0
     ymax        = y1min + (y1max - y1min)/2.0
     zmax        = z1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+5:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #6th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+6:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #7th Children
     xmin        = x1min
     xmax        = x1min + (x1max - x1min)/2.0
     ymin        = y1min + (y1max - y1min)/2.0
     ymax        = y1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+7:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     #8th Children
     xmin        = x1min + (x1max - x1min)/2.0
     xmax        = x1max
     (n,xc,yc,zc,xindex)= mass_center_and_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index)
     mass_tree.update({path+8:(xmin,xmax,ymin,ymax,zmin,zmax,n,xc,yc,zc,xindex)})
     return mass_tree

def add_index(path):
    path = path*10 + 1
    return path

def continue_index(path):
    if (path%10 < 8):
       path  = path+1
    elif(path == 8):
       path  = 0 
    else:
       path  = path/10
       path  = continue_index(path)
    return path

def build_mass_tree(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,n):
    #Loading Data
    (x , y , z)     = put_in_order_data(x,y,z)
    print "data reading done"
    #Initializing Data
    mass_tree       = ini_tree(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z)
#    mass_num        = 0
    print "Initializing done"
    path            = 1
    while(path != 0):
       n            = mass_tree[path][6]
       if (n > 1):
          mass_tree = creat_mass_children(mass_tree,path,x,y,z)
          path      = add_index(path)
#       elif (n == 1):
#          path      = continue_index(path)
#          mass_num += 1
#          print mass_num
       else:
          path      = continue_index(path)
    return (mass_tree , x , y , z)

###
### We are not going to use this function in this implimentation!
###
def find_path(mass_tree,x,y,z):
    path = 0
    for path, v in mass_tree.iteritems():
       if (v[7] == x):
          if (v[8] == y):
             if (v[9] == z):
                 return path
                 break

def force_calculator_OCT(mass_tree,x,y,z,i):
    G      = 6.67e-11 # Gravitational constant
    mass   = 1.0      # Mass of each particle
    angle  = 1.0      # Critical Angle
    #Initializing
    Fx     = 0.0
    Fy     = 0.0
    Fz     = 0.0
    path   = 1
    while(path != 0):
       mass_num= mass_tree[path][6]
       xc      = mass_tree[path][7]
       yc      = mass_tree[path][8]
       zc      = mass_tree[path][9]
       r       = sqrt((xc - x[i])**2 + (yc - y[i])**2 + (zc - z[i])**2) + 10e-20
#       if ((r > 0.0) and (mass_num != 0)):
       d       = mass_tree[path][1] - mass_tree[path][0]
       theta   = d/r
       if ((theta < angle) or (mass_num < 2)): #(mass_num == 1)):
          Fx  += G*mass*mass*mass_num*(xc - x[i])/(r**3)
          Fy  += G*mass*mass*mass_num*(yc - y[i])/(r**3)
          Fz  += G*mass*mass*mass_num*(zc - z[i])/(r**3)
          path = continue_index(path)
       else:
          path = add_index(path)   
#       else:
#          path    = continue_index(path)
    return (Fx,Fy,Fz)

