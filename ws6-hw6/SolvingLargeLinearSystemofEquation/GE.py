# this requires numpy get it from http://numpy.sf.net

from copy import deepcopy
from numpy import *

# this function, swapRows, was adapted from
# Numerical Methods Engineering with Python, Jean Kiusalaas

def swapRows(v,i,j):
	"""Swaps rows i and j of vector or matrix [v]."""
	if len(v) == 1:
		v[i],v[j] = v[j],v[i]
	else:
		temp = v[i].copy()
		v[i] = v[j]
		v[j] = temp
	
def pivoting(a, b):
	"""changes matrix A by pivoting"""
	
	n = len(b)
	
	for k in range(0, n-1):
		p = int(argmax(abs(a[k:n, k]))) + k
		if (p != k):
			swapRows(b, k, p)
			swapRows(a,k,p)

def gauss(a, b, t=1.0e-9, verbose=False):
	""" Solves [a|b] by gauss elimination"""
	
	n = len(b)
	
	# make copies of a and b so as not to change the values in the arguments
	tempa = deepcopy(a)
	tempb = deepcopy(b)
	
	# check if matrix is singular
	if abs(linalg.det(tempa)) < t:
		print "asn"
		return -1
	
	pivoting(tempa, tempb)

	for k in range(0,n-1):	
		for i in range(k+1, n):
			if tempa[i,k] != 0.0:
				m = tempa[i,k]/tempa[k,k]
				if verbose:
				    print "m =", m
				tempa[i,k+1:n] = tempa[i,k+1:n] - m * tempa[k,k+1:n]
				tempb[i] = tempb[i] - m * tempb[k]
	
	# Back substitution
	for k in range(n-1,-1,-1):
		tempb[k] = (tempb[k] - dot(tempa[k,k+1:n], tempb[k+1:n]))/tempa[k,k]

	return tempb

def residue(a, b, c):
    """Calculates the residue of a system solved by gauss elimination"""
    n = len(b)

    t = a * c # t is the A with the values of x replaced (an [n x n] matrix) 

    s = []
    for i in range(0, n):
        s.append(sum(t[i])) # s is the solution


    res = b - s # res is the residue

    return res

# 1) dA = Somar uma pequena quantidade a A
# 2) dB = Somar uma pequena quantidade a B
# 3) A.dX = dB - dA.x0

def disturb(a, b, ai, bi):
    """Calculates the extarnal disturbance cause by da and db"""
    
    x = gauss(a, b)
    
    n = len(b)
    
    da = array([1])
    da.resize(n, n)
    da.fill(ai)
    
    print "da = ", da
    
    db = array([1])
    db.resize(1, n)
    db.fill(bi)
    
    print "db = ", db
    
    t = db - da*x
    nx = gauss(a, t)

    return nx
    


#a = array([[1.0, 2.0, 0.0],[-1.0, 2.0, 3.0],[1.0, 4.0, 1.0]])
#b = array([3.0, -1.0, 4.0])
#a = array([[-1.414214, 2, 0],[1, -1.414214, 1], [0, 2, -1.414214]])
#b = array([1.0,1.0,1.0])
#a = array([[2.0, 2.0, 2.0],[1.0, 1.0, 5.0], [2.0, 5.0, 1.0]])
#b = array([6.0, 7.0, 8.0])
#a = array([[1.001, 2.001, 3.001],[0.999, 2.0, 2.999], [1.002, 1.999, 2.999]])
#b = array([4.003, 4.001, 3.999])

#a = array([[0.7, 8.0, 3.0],[-6.0, 0.45, -0.25], [8.0, -3.1, 1.05]])
#b = array([2.3, 3.5, 10.3])

#x = gauss(a, b)
#print "Solution = ", x

#sol = linalg.solve(a, b)
#print "linalg Solution = ", sol

#y = residue(a, b, x)
#print "Residue = ", y

#u = gauss(a, y)
#print "Residue destribution = ", u

#z = gauss(a, b+y)
#print "New Solution (with added residue) = ", z

#y2 = residue(a, b+y, z)
#print "Residule of new solution = ", y2


#if linalg.norm(y2) < linalg.norm(y):
#    print "New solution has a smaller residue."
#else:
#    print "Original solution has a smaller residue."

#nx = disturb(a, b, 0.5, 0.5)
#print "Disturbed Solution = ", nx

