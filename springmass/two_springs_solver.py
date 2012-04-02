#! /usr/bin/env python

#
# two_springs_solver.py
#
"""
Use ODEINT to solve the differential equations defined by the vector field
in two_springs.py.
"""

import sys
import numpy
from scipy.integrate import odeint
import two_springs

# Parameter values, this could be enter by user as a parameter
# Masses:
m1 = 1.0
m2 = 10.0
# Spring constants
k1 = 10.0
k2 = 1.0
# Natural lengths
L1 = 1.0
L2 = 1.0
# Friction coefficients
b1 = 0.0
b2 = 0.0

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
x1 = 0.0
y1 = 0.0
x2 = 10.0
y2 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 2500

# Set output file using user parameter
if len(sys.argv) > 1:
	f = open(sys.argv[1],'w')
else:
	f = sys.stdout
	
# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
# TODO: How this sentence work?
#t = [stoptime*float(i)/(numpoints-1) for i in range(numpoints)]

# linspace is easier to understand
t = numpy.linspace(0, stoptime, numpoints) # time

# Pack up the parameters and initial conditions:
p = [m1,m2,k1,k2,L1,L2,b1,b2]
w0 = [x1,y1,x2,y2]

# Call the ODE solver, he needs:
# the vector field, 
# the initial condition,
# the evaluate points,
# and arguments for vector field, right?
wsol = odeint(two_springs.vectorfield,w0,t,args=(p,),atol=abserr,rtol=relerr)

# Print the solution,
for t1,w1 in zip(t,wsol):
	f.write('{0} {1} {2} {3} {4}\n'.format(t1,w1[0],w1[1],w1[2],w1[3]))
