import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import time
import os

os.chdir('/home/antony/mycodes/li-rinzel/data/anup')

run = time.strftime("%Y%m%d-%H%M%S")


def PMCA(z, t):
    m0 = z[0]
    m1 = z[1]
    m2 = z[2]
    Ca_cyt = 0.035
    kf1 = 15
    kb1 = 20
    kf2 = 20
    kf3 = 100
    k1 = 0.6
    PMCA_density = 400
    J_pmca = (-m0*Ca_cyt*kf1) + (m1*kb1) + (m0*k1)
    dm0dt = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m2 * kf3)
    dm1dt = (-m1 * kb1) + (m0 * Ca_cyt * kf1) + (-m1 * kf2)
    dm2dt = (-m2 * kf3) + (m1 * kf2)
    dzdt = [dm0dt, dm1dt, dm2dt]
    return dzdt

def cyt_buffer(z,t):
    c0 = z[0]
    c1 = z[1]
    Ca_cyt = 0.035
    kf = 60
    kb = 1200
    buffer_conc = 50
    J_cytB = (-Ca_cyt*c0*kf) + c1*kb
    dc0dt = (-c0*Ca_cyt*kf) + (c1*kb)
    dc1dt = (-c1*kb) + (c0*Ca_cyt*kf)
    dzdt = [dc0dt, dc1dt]
    return dzdt

def SERCA(Ca_cyt):
    V_max = 250
    ka = 100
    n = 2
    J_serca = V_max * ((Ca_cyt**2)/((Ca_cyt **2)+(ka**2)))

def ER_leak(Ca_ER, Ca_Cyto):
    k_leak = 0.2
    J_ERleak = k_leak*(Ca_ER - Ca_Cyto)
    return J_ERleak
