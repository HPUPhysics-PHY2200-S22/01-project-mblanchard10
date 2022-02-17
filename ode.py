# Program: ODE solvers library (ode.py)
# J Wang, Computational modeling and visualization with Python
#

import numpy as np      # numpy

def Euler(diffeq, y0, t, h): # uses docstring """..."""
    """ Euler's method for n ODEs:
        Given y0 at t, returns y1 at t+h """
    dydt = diffeq(y0, t)     # get {dy/dt} at t
    return y0 + h*dydt       # Euler method on a vector  @\lbl{line:vecop}@
 
def RK2(diffeq, y0, t, h):
    """ RK2 method for ODEs:
        Given y0 at t, returns y1 at t+h """
    k1 = h*diffeq(y0, t)                # get dy/dt at t first
    k2 = h*diffeq(y0+0.5*k1, t + h/2)   # get dy/dt at t+h/2,
    return y0 + k2                      # calc. y1 = y(t+h)
    
def RK4(diffeq, y0, t, h):
    """ RK4 method for ODEs:
        Given y0 at t, returns y1 at t+h """
    k1 = h*diffeq(y0, t)                    # dy/dt at t
    k2 = h*diffeq(y0+0.5*k1, t + h/2.)      # dy/dt at t+h/2
    k3 = h*diffeq(y0+0.5*k2, t + h/2.)      # dy/dt at t+h/2
    k4 = h*diffeq(y0+k3, t + h)             # dy/dt at t+h
    return y0 + (k1+k4)/6.0 + (k2+k3)/3.0
