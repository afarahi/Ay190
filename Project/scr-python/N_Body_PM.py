import sys
from pylab         import *
from scipy         import *
from numpy         import *
from numpy.fft     import fftn
from math          import sqrt, pi, sin, cos
from cmath         import exp

def mass_counter(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,index0=0):
    n     = 0
    for i in range(len(x)):
       if ((x[i]>=xmin) and (x[i]<xmax)):
          if((y[i]>=ymin) and (y[i]<ymax)):
             if((z[i]>=zmin) and (z[i]<zmax)):
                n += 1
    return n

def density_kernel_I(x0,y0,z0,x,y,z):
    mass    = 1.0      # Mass of each particle
    density = 0.0
    q       = 0.1
    for i in range(len(x)):
       R = sqrt((x[i]-x0)**2 + (y[i]-y0)**2 + (z[i]-z0)**2)
       if (R < q):
         if (R < 0.5*q):
            density += 8.0/pi * (1.0 - 6.0*R**2 + 0.6*R**3)
         else:
            density += 16.0/pi * (1.0 - R)**3
    return density*mass

def density_CIC(x0,y0,z0,x,y,z,dx):
    mass    = 1.0      # Mass of each particle
    density = 0.0
    for i in range(len(x)):
       if ((abs(x[i] - x0)) < dx and (abs(y[i] - y0) < dx) and (abs(z[i] - z0) < dx) ):
            density += mass*((1.0 - abs(x[i] - x0)/dx)*(1.0 - abs(y[i] - y0)/dx)*(1.0 - abs(z[i] - z0)/dx))
    return density

def mesh_the_volume(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,n):
    mass     = 1.0      # Mass of each particle
    dx       = (xmax - xmin)/n
    dy       = (ymax - ymin)/n
    dz       = (zmax - zmin)/n
    volume   = dx*dy*dz
    R        = zeros([n,n,n,4])
    l        = 0
    for i in range(n):
       for j in range(n):
          for k in range(n):
             x0min   = xmin + i*dx
             y0min   = ymin + j*dy
             z0min   = zmin + k*dz
             x0max   = xmin + (i+1)*dx
             y0max   = ymin + (j+1)*dy
             z0max   = zmin + (k+1)*dz
             R[i,j,k,1]    = (x0max + x0min)/2.0
             R[i,j,k,2]    = (y0max + y0min)/2.0
             R[i,j,k,3]    = (z0max + z0min)/2.0
#             R[i,j,k,0]    = mass*mass_counter(x0min,x0max,y0min,y0max,z0min,z0max,x,y,z)/volume #NGP
             R[i,j,k,0]    = density_CIC(R[i,j,k,1],R[i,j,k,2],R[i,j,k,3],x,y,z,dx)/volume #CIC
#             R[i,j,k,0]    = density_kernel_I(R[i,j,k,1],R[i,j,k,2],R[i,j,k,3],x,y,z)/volume #Kernel
             l += 1
    return R

def trilinear_interp(input_array, indices):
    """Evaluate the input_array data at the indices given"""
    output = np.empty(indices[0].shape)
    x_indices = indices[0]
    y_indices = indices[1]
    z_indices = indices[2]
    for i in np.ndindex(x_indices.shape):
        x0 = np.floor(x_indices[i])
        y0 = np.floor(y_indices[i])
        z0 = np.floor(z_indices[i])
        x1 = x0 + 1
        y1 = y0 + 1
        z1 = z0 + 1
        #Check if xyz1 is beyond array boundary:
        if x1 == input_array.shape[0]:
            x1 = x0
        if y1 == input_array.shape[1]:
            y1 = y0
        if z1 == input_array.shape[2]:
            z1 = z0
        x = x_indices[i] - x0
        y = y_indices[i] - y0
        z = z_indices[i] - z0
        output[i] = (input_array[x0,y0,z0]*(1-x)*(1-y)*(1-z) +
                 input_array[x1,y0,z0]*x*(1-y)*(1-z) +
                 input_array[x0,y1,z0]*(1-x)*y*(1-z) +
                 input_array[x0,y0,z1]*(1-x)*(1-y)*z +
                 input_array[x1,y0,z1]*x*(1-y)*z +
                 input_array[x0,y1,z1]*(1-x)*y*z +
                 input_array[x1,y1,z0]*x*y*(1-z) +
                 input_array[x1,y1,z1]*x*y*z)

    return output

def potential_field(k,density):
    G        = 6.67e-11 # Gravitational constant
    p        = - 4.0 * pi * G * density / k
    return p

def k_space(n,L):
    k = zeros(n)
    if (mod(n,2) == 0):
      for i in range(n/2):
         k[i]     = i * (2*pi/L)
         k[n-1-i] = (-i-1) * (2*pi/L)
      k[n/2]  = 0.0
    else:
      print "For odd numbers it is not working ... ! "
      sys.exit()
    return k

def force_calculator_MP(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,n):
   G         = 6.67e-11 # Gravitational constant
   n3        = n*n*n
   L         = xmax - xmin
   #Meshing
   R         = mesh_the_volume(xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,n)
   #Fourier Space
   density_k = fftn(R[:,:,:,0])
   #Density in real space (just for checking the fourier transformation)
   #density_x = ifftn(density_k)
   #print density_x.round()
   #Force in k-space
   k_s       = k_space(n,L)
   Fx_k      = zeros([n,n,n],complex)
   for i in range(n):
      for j in range(n):
         for k in range(n):
            k2               = k_s[i]**2 + k_s[j]**2 + k_s[k]**2
            if (k2 != 0.0):
               Fx_k[i,j,k]   = 4.0j*k_s[i]*pi*G*density_k[i,j,k]/k2 
   F             = ifftn(Fx_k)
   ### Interpolation
   #####
   Force     = zeros([len(x)],complex)
   dx        = (xmax - xmin)/n
   for l in range(len(x)):
      mm = 0.0
      for i in range(n):
         for j in range(n):
            for k in range(n):
#               DISTANCE = sqrt((R[i,j,k,1]-x[l])**2 + (R[i,j,k,2]-y[l])**2 + (R[i,j,k,3]-z[l])**2)
#               if (DISTANCE <= dx):
#                 Force[l] += F[i,j,k]*DISTANCE
#                 mm       += DISTANCE
#      Force[l] = Force[l]/mm
               if (abs(R[i,j,k,1]-x[l]) <= dx/2.0 and abs(R[i,j,k,2]-y[l]) <= dx/2.0 and abs(R[i,j,k,3]-z[l]) <= dx/2.0):
                   Force[l] = F[i,j,k]
   return Force
