#! /usr/bin/env python

import numpy
import pylab as pl
from scipy.integrate import odeint

# Definition of parameters
abserr = 1.0e-8
relerr = 1.0e-6
a = 1.
b = 0.1
c = 1.5
d = 0.75
stoptime = 10
numpoints = 1000

prey_0, predator_0 = 10, 10

def dX_dt(X, t=0):
	""" 
	Return the growth rate of fox and rabbit populations. 
	"""
	f = [a*X[0] - b*X[0]*X[1], 
		-c*X[1] + d*b*X[0]*X[1]]
	return f

t = numpy.linspace(0, stoptime, numpoints) # time
X0 = [prey_0, predator_0]            # initials conditions

X = odeint(dX_dt, X0, t, atol=abserr, rtol=relerr)

# We could also transpose X ...
# prey, predator = X.T

# But this is more clear to python newbies (as me)
prey = X[:,0]
predator = X[:,1]
f1 = pl.figure()

pl.plot(t, prey, 'r-', label='Prey')
pl.plot(t, predator  , 'b-', label='Predator')

pl.grid()
pl.legend(loc='best')
pl.xlabel('time')
pl.ylabel('population')
pl.title('Evolution of prey and predator populations')

# we construct values for generate 
# multiple initial values, 
# see X0 below
values = numpy.linspace(0.3, 0.9, 5) 

# equilibrium points
X_f0 = [     0. ,  0.]
X_f1 = [ c/(d*b), a/b]

# colors for each trajectory
vcolors = pl.cm.autumn_r(numpy.linspace(0.3, 1., len(values)))  

f2 = pl.figure()

# plot trajectories
for v, col in zip(values, vcolors):
	# starting point
	X0 = v * numpy.array(X_f1)
	X = odeint(dX_dt, X0, t)
	pl.plot( X[:,0], X[:,1], lw=3.5*v, color=col, label='X0=(%.f, %.f)' % ( X0[0], X0[1]))

# define a grid and compute direction at each point
ymax = pl.ylim(ymin=0)[1]                        # get axis limits
xmax = pl.xlim(xmin=0)[1]
nb_points = 20

x = numpy.linspace(0, xmax, nb_points)
y = numpy.linspace(0, ymax, nb_points)

X1 , Y1  = numpy.meshgrid(x, y)                       # create a grid
DX1, DY1 = dX_dt([X1, Y1])                      # compute growth rate on the gridt

M = (numpy.hypot(DX1, DY1))                           # Norm of the growth rate 
M[ M == 0] = 1.                                 # Avoid zero division errors 
DX1 /= M                                        # Normalize each arrows
DY1 /= M

# Draw direction fields, using matplotlib 's quiver function
# I choose to plot normalized arrows and to use colors to give information on
# the growth speed
pl.title('Trajectories and direction fields')
Q = pl.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=pl.cm.jet)
pl.xlabel('Prey')
pl.ylabel('Predator')
pl.legend()
pl.grid()
pl.xlim(0, xmax)
pl.ylim(0, ymax)

def IF(X):
	u, v = X
	return u**(c/a) * v * numpy.exp( -(b/a)*(d*u+v) )

# Plot iso contours
nb_points = 80                              # grid size
x = numpy.linspace(0, xmax, nb_points)
y = numpy.linspace(0, ymax, nb_points)
X2 , Y2  = numpy.meshgrid(x, y)                   # create the grid
Z2 = IF([X2, Y2])                           # compute IF on each point
f3 = pl.figure()
CS = pl.contourf(X2, Y2, Z2, cmap=pl.cm.Purples_r, alpha=0.5)
CS2 = pl.contour(X2, Y2, Z2, colors='black', linewidths=2. )
pl.clabel(CS2, inline=1, fontsize=16, fmt='%.f')
pl.grid()
pl.xlabel('Prey')
pl.ylabel('Predator')
pl.ylim(1, ymax)
pl.xlim(1, xmax)
pl.title('IF contours')

pl.show()


