from pylab import *
import sys,math

def analytic(t,t0,x,x0,D):
    y = (t0/(t0+t))**(1.0/2.0)*exp(-(x - x0)**2 / (4.0 * D * (t + t0)) )       
#    print "Implement me!!!"
    return y

n   = 200
x   = linspace(0,100,n)
dx  = x[1]-x[0]
dx2 = dx*dx

x0  = 50.0
D   = 1.0
t0  = 1.0
t   = 0.0
nt  = 3000


# initial conditions
y = exp( -(x-x0)**2 / (4.0*D*t0) )

dt = 1.0 * dx**2 / (2.0 * D) #"SET ME"

def calc_rhs(D,y,dx):
    RHS      = zeros(len(y))
    for i in range(1,len(y)-1):
       RHS[i]=  D*(y[i+1]-2.0*y[i]+y[i-1])/(dx*dx) 
#    print "Implement me!"
    return RHS

def set_bound(y):
    y[0] = y[1]
    y[n-1] = y[n-2]
    return y

ion()
plot(x,y,"r-")
plot(x,analytic(t,t0,x,x0,D),"bx-")
show()

for it in range(nt):
    RHS = calc_rhs(D,y,dx)
#    y[1:-1] = y[1:-1] + dt * RHS[1:-1] #"Implement me"
    for i in range(1,len(y)-1):
      y[i] = y[i] + dt*RHS[i]

    y = set_bound(y)
    t = t + dt

    if it % 10 == 0:
        print it
        clf()
        plot( (0,0), (0,1.0), color="white")
        plot(x,y,'rx-')
        plot(x,analytic(t,t0,x,x0,D),"bx-")
        draw()

ioff()
plot(x,y,"r-")
show()


