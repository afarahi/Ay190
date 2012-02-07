import sys
from pylab import *
from math  import exp, sqrt, fmod
from scipy import *

def apply_bcs(x):
    # apply boundary conditions
    # you need to fill in code
    sigma = sqrt(15.0)
    x0    = 30.0
    y     = exp(-(x-x0)**2/(2.0*sigma**2))  
    return y

def analytic(x,t):
    sigma = sqrt(15.0)
    x0    = 30.0
    y     = exp(-(x-x0-v*t)**2/(2.0*sigma**2))  
    return y

def central_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew[i] = yold[i] - (v*dt/(2.0*dx))*(yold[i+1] - yold[i-1])
    ynew[0]   = y[1]
    ynew[n-1] = y[n-2]
    return ynew

def upwind_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew[i] = yold[i] - (v*dt/(dx))*(yold[i+1] - yold[i])
    ynew[0]   = y[1]
    ynew[n-1] = y[n-2]
    return ynew

def downwind_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew[i] = yold[i] - (v*dt/(dx))*(yold[i] - yold[i-1])
    ynew[0]   = y[1]
    ynew[n-1] = y[n-2]
    return ynew

def Lax_Friedrichs_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew[i] = (yold[i+1]+yold[i-1])/2.0 - (v*dt/(2.0*dx))*(yold[i+1] - yold[i-1])
    ynew[0]   = y[1]
    ynew[n-1] = y[n-2]
    return ynew

def MacCormack_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew1= zeros(n)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew1[i] = yold[i] - (v*dt/(dx))*(yold[i+1] - yold[i])    
    ynew1[0]   = y[1]
    ynew1[n-1] = y[n-2]
    for i in range(1,n-1):
        ynew[i] = (ynew1[i]+yold[i])/2.0 - (v*dt/(2.0*dx))*(ynew1[i] - ynew1[i-1])
    ynew[0]   = ynew1[1]
    ynew[n-1] = ynew1[n-2]
    return ynew

def Lax_Wendroff_sol(yold,v,dx,dt):
    n    = len(yold)
    ynew1= zeros(n)
    ynew = zeros(n)
    for i in range(1,n-1):
        ynew1[i] = (yold[i+1]+yold[i])/2.0 + (v*dt/(2.0*dx))*(yold[i+1] - yold[i])    
    ynew1[0]   = y[1]
    ynew1[n-1] = y[n-2]
    for i in range(1,n-1):
        ynew[i] = yold[i] + (v*dt/dx)*(ynew1[i+1] - ynew1[i])
    ynew[0]   = ynew1[1]
    ynew[n-1] = ynew1[n-2]
    return ynew


# parameters
n = 1000
v = 0.1
# set up the grid here. Use a decent number of zones;
# perhaps to get a dx of 0.1
x     = linspace(0.0,100.0,n+1)
dx    = x[1] - x[0]
y     = zeros(n+1)
cfl   = 1.0
dt    = 1.0
t     = 0.0
v     = 0.1


#set up initial conditions
y = apply_bcs(x)

# evolve (and show evolution)
figure()
plot(x,y,'x-') # numerical data
plot(x,analytic(x,t),'r-') # analytic data
show()

yold  = y
ntmax = 1000

'''for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = upwind_sol(yold,v,dx,dt)
    yan  = analytic(x,t)  
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    semilogy(t,abs(err),'r.')
semilogy(t,abs(err),'r.', label='Upwind scheme')

#    clf()
#    plot(x,y,'b-')
#    plot(x,yan,'r-')
#    draw()

t     = 0.0

#set up initial conditions
y = apply_bcs(x)

for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = downwind_sol(yold,v,dx,dt)
    yan  = analytic(x,t)  
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    plot(t,abs(err),'b.')
plot(t,abs(err),'b.', label='Downwind scheme')

t     = 0.0

#set up initial conditions
y = apply_bcs(x)

for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = central_sol(yold,v,dx,dt)
    yan  = analytic(x,t)
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    semilogy(t,abs(err),'m.')
semilogy(t,abs(err),'m.', label='Central scheme')

t     = 0.0

#set up initial conditions
y = apply_bcs(x)

for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = Lax_Friedrichs_sol(yold,v,dx,dt)
    yan  = analytic(x,t)
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    semilogy(t,abs(err),'r.')
semilogy(t,abs(err),'r.', label='Lax_Friedrichs Method')

t     = 0.0

#set up initial conditions
y = apply_bcs(x)

for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = MacCormack_sol(yold,v,dx,dt)
    yan  = analytic(x,t)
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    semilogy(t,abs(err),'b.')
semilogy(t,abs(err),'b.', label='MacCormack Method')


t     = 0.0

#set up initial conditions
y = apply_bcs(x)

for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = Lax_Wendroff_sol(yold,v,dx,dt)
    yan  = analytic(x,t)
    err  = (max(y)-max(yan))/max(yan)
    figure(1)
    semilogy(t,abs(err),'g.')
semilogy(t,abs(err),'g.', label='Lax-Wendroff Method')


xlabel('time')
xlabel('relative error')
legend(loc=2)
savefig('plot5-2.pdf')
'''

ion()
figure()
for it in range(ntmax):
    print it
    t    = t + dt
    # save previous and previous previous data
    yold = y
    # get new data; ideally just call a function
    y    = MacCormack_sol(yold,v,dx,dt)
    if (fmod(it,30) == 0):
       yan  = analytic(x,t)
       clf()
       plot(x,y,'b-')
       plot(x,yan,'r-')
       draw()


