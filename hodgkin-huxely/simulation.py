# This simulations uses the parameters given my Hodgkin-Huxely corrected to a membrane potential of -65 mv
import numpy as np

# Value of the membrane Capacitance

Cm = 1

# Value for reversal potentials

ENa = 50
EK = -77
EL = -54.4

# Value for conductance

gNa = 120
gK = 36
gL = 0.3

u = -65


# Equations for the gating variables

alpha_n = (0.1 - 0.01*u) / (np.exp(1-0.1*u) - 1)

beta_n = 0.125*np.exp(-u/80)

alpha_m = (2.5 - 0.1*u) / (np.exp(2.5-0.1*u) - 1) 

beta_m = 4*np.exp(-u/18)

alpha_h =  0.07*np.exp(-u/20)

beta_h =  1 / (np.exp(3-0.1*u) + 1)

# Compute values for n,m,h

def Iin():
    return 1

def model(alpha_m, beta_m,alpha_n, beta_n,alpha_h, beta_h):
    dydt = [alpha_m*(1-m)]
    return dydt

