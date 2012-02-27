# hardcoded for 3 inner boundary points
def apply_bcs_spherical(hyd):
    hyd.rho[0]   = hyd.rho[5]
    hyd.rho[1]   = hyd.rho[4]
    hyd.rho[2]   = hyd.rho[3]
    hyd.vel[0]   = -hyd.vel[5]
    hyd.vel[1]   = -hyd.vel[4]
    hyd.vel[2]   = -hyd.vel[3]
    hyd.eps[0]   = hyd.eps[5]
    hyd.eps[1]   = hyd.eps[4]
    hyd.eps[2]   = hyd.eps[3]
    hyd.press[0] = hyd.press[5]
    hyd.press[1] = hyd.press[4]
    hyd.press[2] = hyd.press[3]
    
    hyd.rho[hyd.n-hyd.g:hyd.n-1]   = hyd.rho[hyd.n-hyd.g-1]
    hyd.vel[hyd.n-hyd.g:hyd.n-1]   = hyd.vel[hyd.n-hyd.g-1]
    hyd.eps[hyd.n-hyd.g:hyd.n-1]   = hyd.eps[hyd.n-hyd.g-1]
    hyd.press[hyd.n-hyd.g:hyd.n-1] = hyd.press[hyd.n-hyd.g-1]
    return hyd
