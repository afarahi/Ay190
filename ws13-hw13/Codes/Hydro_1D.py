#!/usr/bin/env python
import sys,math
from pylab import *

# basic parameters
gamma = 5.0/3.0
cfl   = 0.5
dt    = 1.0e-5
dtp   = dt
reconstruction_type = 'pc' # minmod, mc, pc
# use nzones
nzones= 207
# run until time 0.2
tend  = 0.2

##################################### class definition
class mydata:
    def __init__(self,nzones):
        self.x      = zeros(nzones) # cell centers
        self.xi     = zeros(nzones) # cell LEFT interfaces
        self.rho    = zeros(nzones)
        self.rhop   = zeros(nzones)
        self.rhom   = zeros(nzones)
        self.vel    = zeros(nzones)
        self.velp   = zeros(nzones)
        self.velm   = zeros(nzones)
        self.eps    = zeros(nzones)
        self.epsp   = zeros(nzones)
        self.epsm   = zeros(nzones)
        self.press  = zeros(nzones)
        self.pressp = zeros(nzones)
        self.pressm = zeros(nzones)
        self.q      = zeros((3,nzones)) # conserved quantities
        self.qp     = zeros((3,nzones))  
        self.qm     = zeros((3,nzones))  
        self.n      = nzones
        self.g      = 3

    def setup_grid(self,xmin,xmax):
        dx   = (xmax - xmin) / (self.n - self.g*2 - 1)
        xmin = xmin - self.g*dx
        xmax = xmax + self.g*dx
        for i in range(self.n):
            # cell centers
            self.x[i]  = xmin + (i)*dx
        # cell LEFT interfaces
        for i in range(self.n):
            self.xi[i] = self.x[i] - 0.5*dx

    def setup_ID(self):
        # Shocktube initial data
        rchange = (self.x[self.n-self.g-1] - self.x[self.g]) / 2.0
        rho1    = 1.0
        rho2    = 0.1
        press1  = 1.0
        press2  = 0.125
        for i in range(self.n):
            if self.x[i] < rchange:
                self.rho[i]   = rho1
                self.press[i] = press1
                self.eps[i]   = press1 / (rho1) / (gamma - 1.0)
                self.vel[i]   = 0.0
            else:
                self.rho[i]   = rho2
                self.press[i] = press2
                self.eps[i]   = press2 / (rho2) / (gamma - 1.0)
                self.vel[i]   = 0.0

##################################### some basic functions
def prim2con(rho,vel,eps):
    q = zeros((3,len(rho)))
    q[0,:] = rho[:]
    q[1,:] = rho[:] * vel[:]
    q[2,:] = rho[:] * eps[:] + \
                  0.5 * rho[:] * vel[:]**2
    return q

def con2prim(q):
    rho = q[0,:]
    vel = q[1,:] / rho
    eps = q[2,:] / rho - 0.5*vel**2
    press = eos_press(rho,eps,gamma)
    return (rho,eps,press,vel)

############# boundary conditions
def apply_bcs(hyd):
    hyd.rho[0:hyd.g-1]   = hyd.rho[hyd.g]
    hyd.vel[0:hyd.g-1]   = hyd.vel[hyd.g]
    hyd.eps[0:hyd.g-1]   = hyd.eps[hyd.g]
    hyd.press[0:hyd.g-1] = hyd.press[hyd.g]

    hyd.rho[hyd.n-hyd.g:hyd.n-1]   = hyd.rho[hyd.n-hyd.g-1]
    hyd.vel[hyd.n-hyd.g:hyd.n-1]   = hyd.vel[hyd.n-hyd.g-1]
    hyd.eps[hyd.n-hyd.g:hyd.n-1]   = hyd.eps[hyd.n-hyd.g-1]
    hyd.press[hyd.n-hyd.g:hyd.n-1] = hyd.press[hyd.n-hyd.g-1]
    return hyd

############# minmod function 
def minmod(a,b):
    if(a*b < 0):
        mm = 0.0
    elif(abs(a)<abs(b)):
        mm=a
    else:
        mm=b
    return mm

############# monotonized central function 
def monotonized_central(a,b):
    if(a*b < 0):
        mm = 0.0
    else:
        mm= signum(min( 2*abs(a) , 2*abs(b) , 0.5 *(abs(a) + abs(b)) ) , a + b )
    return mm

############# minmod function 
def tvd_minmod_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
         du_up = (f[i] - f[i-1]) / (x[i] - x[i-1])
         du_do = (f[i+1] - f[i]) / (x[i+1] - x[i])
         delta = minmod(du_up,du_do)
         fp[i] = f[i] + delta*(xi[i+1] - x[i])
         fm[i] = f[i] - delta*(x[i]   - xi[i])
#        print "implement me! tvd_minmod_reconstruction :-)"
#        sys.exit()
    return (fp,fm)

############# signum functions
def signum(x,y):
    if(y >= 0):
        return abs(x)
    else:
        return -abs(x)

############# mc reconstruction
def tvd_mc_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
         du_up = (f[i] - f[i-1]) / (x[i] - x[i-1])
         du_do = (f[i+1] - f[i]) / (x[i+1] - x[i])
         delta = monotonized_central(du_up,du_do)
         fp[i] = f[i] + delta*(xi[i+1] - x[i])
         fm[i] = f[i] - delta*(x[i]   - xi[i])
    return (fp,fm)


############# reconstruction top level function
def reconstruct(hyd,type):
    if(type=='pc'):
        # piecewise constant reconstruction 
        for i in range(hyd.g-1,hyd.n-hyd.g+1):
            hyd.rhop[i] = hyd.rho[i]
            hyd.rhom[i] = hyd.rho[i]
            hyd.epsp[i] = hyd.eps[i]
            hyd.epsm[i] = hyd.eps[i]
            hyd.velp[i] = hyd.vel[i]
            hyd.velm[i] = hyd.vel[i]

    elif(type=='minmod'):
        (hyd.rhop,hyd.rhom) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)

    elif(type=='mc'):
        (hyd.rhop,hyd.rhom) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)
                
    else:
        print "reconstruction type not known; abort!"
        sys.exit()

    hyd.pressp = eos_press(hyd.rhop,hyd.epsp,gamma)
    hyd.pressm = eos_press(hyd.rhom,hyd.epsm,gamma)

    hyd.qp = prim2con(hyd.rhop,hyd.velp,hyd.epsp)
    hyd.qm = prim2con(hyd.rhom,hyd.velm,hyd.epsm)

    return hyd

############# equation of state
def eos_press(rho,eps,gamma):
    press  = (gamma - 1.0) * rho * eps
    return press

def eos_cs2(rho,eps,gamma):
    prs    = (gamma - 1.0) * rho * eps
    dpde   = (gamma - 1.0) * rho
    dpdrho = (gamma - 1.0) * eps
    cs2    = dpdrho + dpde * prs/(rho+1.0e-30)**2
    return cs2

############# time step calculation
def calc_dt(hyd,dtp):
    cs    = sqrt(eos_cs2(hyd.rho,hyd.eps,gamma))
    dtnew = 1.0
    for i in range(hyd.g,hyd.n-hyd.g):
        dtnew = min(dtnew, (hyd.x[i+1]-hyd.x[i]) / \
                    max(abs(hyd.vel[i]+cs[i]), \
                        abs(hyd.vel[i]-cs[i])))
    dtnew = min(cfl*dtnew,1.05*dtp)
    return dtnew

############# HLLE solver
def hlle(hyd):
    fluxdiff = zeros((3,hyd.n))
    # compute eigenvalues
    evl      = zeros((3,hyd.n))
    evr      = zeros((3,hyd.n))
    smin     = zeros(hyd.n)
    smax     = zeros(hyd.n)
    csp      = sqrt(eos_cs2(hyd.rhop,hyd.epsp,gamma))
    csm      = sqrt(eos_cs2(hyd.rhom,hyd.epsm,gamma))
    for i in range(1,hyd.n-2):
        evl[0,i] = hyd.velp[i]   + csp[i]
        evr[0,i] = hyd.velm[i+1] + csm[i+1]
        evl[1,i] = hyd.velp[i]   
        evr[1,i] = hyd.velm[i+1] 
        evl[2,i] = hyd.velp[i]   - csp[i]
        evr[2,i] = hyd.velm[i+1] - csm[i+1]
        smin[i]  = min(evr[0,i],evl[0,i],evr[1,i],evl[1,i],evr[2,i],evl[2,i],0.0)
        smax[i]  = max(evr[0,i],evl[0,i],evr[1,i],evl[1,i],evr[2,i],evl[2,i],0.0)

    # set up flux left L and right R of the interface
    # at i+1/2
    fluxl = zeros((3,hyd.n))
    fluxr = zeros((3,hyd.n))
    
    for i in range(1,hyd.n-2):
        fluxl[0,i] = hyd.rhop[i]*hyd.velp[i]
        fluxr[0,i] = hyd.rhom[i+1]*hyd.velm[i+1]
        fluxl[1,i] = hyd.rhop[i]*hyd.velp[i]*hyd.velp[i]       + hyd.pressp[i]
        fluxr[1,i] = hyd.rhom[i+1]*hyd.velm[i+1]*hyd.velm[i+1] + hyd.pressm[i+1]
        fluxl[2,i] = (hyd.rhop[i]*hyd.epsp[i]     + 0.5*hyd.rhop[i]*hyd.velp[i]**2     + hyd.pressp[i]  )*hyd.velp[i] 
        fluxr[2,i] = (hyd.rhom[i+1]*hyd.epsm[i+1] + 0.5*hyd.rhom[i+1]*hyd.velm[i+1]**2 + hyd.pressm[i+1])*hyd.velm[i+1]
    
    # solve the Riemann problem for the i+1/2 interface
    ds   = smax - smin
    flux = zeros((3,hyd.n))
    for i in range(hyd.g-1,hyd.n-hyd.g+1):
        flux[0,i] = (smax[i]*fluxl[0,i] - smin[i]*fluxr[0,i] + smin[i]*smax[i]*(hyd.qm[0,i+1]-hyd.qp[0,i])) / ds[i]
        flux[1,i] = (smax[i]*fluxl[1,i] - smin[i]*fluxr[1,i] + smin[i]*smax[i]*(hyd.qm[1,i+1]-hyd.qp[1,i])) / ds[i]
        flux[2,i] = (smax[i]*fluxl[2,i] - smin[i]*fluxr[2,i] + smin[i]*smax[i]*(hyd.qm[2,i+1]-hyd.qp[2,i])) / ds[i]

    # flux differences
    for i in range(hyd.g,hyd.n-hyd.g):
        rm = hyd.xi[i]
        rp = hyd.xi[i+1]
        dxi = 1.0/(rp - rm)
        fluxdiff[0,i] = dxi*( flux[0,i] - flux[0,i-1] )
        fluxdiff[1,i] = dxi*( flux[1,i] - flux[1,i-1] )
        fluxdiff[2,i] = dxi*( flux[2,i] - flux[2,i-1] )

    return fluxdiff

############# RHS calculation
def calc_rhs(hyd,reconstruction_type = 'pc'):
    # reconstruction and prim2con
    hyd      = reconstruct(hyd,reconstruction_type)
    # compute flux differences
    fluxdiff = hlle(hyd)
    # return RHS = - fluxdiff
    return -fluxdiff

########################################################################
# Main program
########################################################################
method = 'pc'

# initialize
hyd = mydata(nzones)

# set up grid
hyd.setup_grid(0.0,1.0)

# set up initial data
hyd.setup_ID()

# get initial timestep
dt = calc_dt(hyd,dt)

# initial prim2con
hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

t = 0.0
i = 0

# display stuff
#ion()
#figure()
#plot(hyd.x,hyd.rho,"r-")
#show()

# main integration loop
while(t < tend):
    
   if(i % 10 == 0):
      # output
      print "%5d %15.6E %15.6E" % (i,t,dt)
#      clf()
#      plot(hyd.x,hyd.rho,"r-")
#      draw()

   # calculate new timestep
   dt     = calc_dt(hyd,dt)

   # save old state
   hydold = hyd
   qold   = hyd.q

   # calc rhs
   k1      = calc_rhs(hyd,method)

   # calculate intermediate step
   hyd.q = qold + 1.0/2.0 * dt * k1
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # boundaries
   hyd = apply_bcs(hyd)

   #calc rhs
   k2 = calc_rhs(hyd,method)
   #apply update
   hyd.q = qold + dt * (0.5 * k1 + 0.5 * k2)
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # apply bcs
   hyd = apply_bcs(hyd)

   # update time
   t = t + dt
   i = i + 1

# Loading Data
A = loadtxt("data.txt")
x = list(A[:,1])
de= list(A[:,2])
# display final result
# ioff()
figure(1)
plot(hyd.x,hyd.rho,"r-", label= 'PC')
plot(x,de,"b-", label = 'Exact')
xlabel('X')
ylabel('Density')

ERR = zeros(hyd.n-7)
for i in range(len(ERR)):
   ERR[i] = abs((hyd.rho[i+3]-de[i])/de[i]) 

figure(2)
plot(x,ERR,"r-",label = 'PC')
xlabel('X')
ylabel('Relative Error')

method = 'minmod'

# initialize
hyd = mydata(nzones)

# set up grid
hyd.setup_grid(0.0,1.0)

# set up initial data
hyd.setup_ID()

# get initial timestep
dt = calc_dt(hyd,dt)

# initial prim2con
hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

t = 0.0
i = 0

# main integration loop
while(t < tend):
    
   if(i % 10 == 0):
      # output
      print "%5d %15.6E %15.6E" % (i,t,dt)
#      clf()
#      plot(hyd.x,hyd.rho,"r-")
#      draw()

   # calculate new timestep
   dt     = calc_dt(hyd,dt)

   # save old state
   hydold = hyd
   qold   = hyd.q

   # calc rhs
   k1      = calc_rhs(hyd,method)

   # calculate intermediate step
   hyd.q = qold + 1.0/2.0 * dt * k1
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # boundaries
   hyd = apply_bcs(hyd)

   #calc rhs
   k2 = calc_rhs(hyd,method)
   #apply update
   hyd.q = qold + dt * (0.5 * k1 + 0.5 * k2)
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # apply bcs
   hyd = apply_bcs(hyd)

   # update time
   t = t + dt
   i = i + 1


ERR = zeros(hyd.n-7)
for i in range(len(ERR)):
   ERR[i] = abs((hyd.rho[i+3]-de[i])/de[i]) 

figure(1)
plot(hyd.x,hyd.rho,"m-",label = 'Minmod')

figure(2)
plot(x,ERR,"m-",label = 'Minmod')




method = 'mc'

# initialize
hyd = mydata(nzones)

# set up grid
hyd.setup_grid(0.0,1.0)

# set up initial data
hyd.setup_ID()

# get initial timestep
dt = calc_dt(hyd,dt)

# initial prim2con
hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

t = 0.0
i = 0

# main integration loop
while(t < tend):
    
   if(i % 10 == 0):
      # output
      print "%5d %15.6E %15.6E" % (i,t,dt)
#      clf()
#      plot(hyd.x,hyd.rho,"r-")
#      draw()

   # calculate new timestep
   dt     = calc_dt(hyd,dt)

   # save old state
   hydold = hyd
   qold   = hyd.q

   # calc rhs
   k1      = calc_rhs(hyd,method)

   # calculate intermediate step
   hyd.q = qold + 1.0/2.0 * dt * k1
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # boundaries
   hyd = apply_bcs(hyd)

   #calc rhs
   k2 = calc_rhs(hyd,method)
   #apply update
   hyd.q = qold + dt * (0.5 * k1 + 0.5 * k2)
   # con2prim
   (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
   # apply bcs
   hyd = apply_bcs(hyd)

   # update time
   t = t + dt
   i = i + 1


ERR = zeros(hyd.n-7)
for i in range(len(ERR)):
   ERR[i] = abs((hyd.rho[i+3]-de[i])/de[i]) 

figure(1)
plot(hyd.x,hyd.rho,"k-",label = 'MC')
legend()
savefig('plot1.pdf')


figure(2)
plot(x,ERR,"k-",label = 'MC')
legend()
savefig('plot2.pdf')

show()
