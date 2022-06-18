import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import math

os.chdir('/home/antony/mycodes/li-rinzel/data/anup')

run = time.strftime("%Y%m%d-%H%M%S")

def model(z,t):
    # Calcium Dynamics
    # ******************************************************************************
    #flux due to PMCA
    m0 = z[]
    m1 = z[]
    m2 = z[]
    Ca_cyt = z[]
    Ca_ER = z[]
    IP3 = z[]
    kf1 = 15
    kb1 = 20
    kf2 = 20
    kf3 = 100
    k1 = 0.6
    dm0dt = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m2 * kf3)
    dm1dt = (-m1 * kb1) + (m0 * Ca_cyt * kf1) + (-m1 * kf2)
    dm2dt = (-m2 * kf3) + (m1 * kf2)
    J_pmca = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m0 * k1)

    # flux from cytoplasmic calcium buffering
    kf = 60
    kb = 1200
    c0 = z[]
    c1 = z[]
    dc0dt = (-c0 * Ca_cyt * kf) + (c1 * kb)
    dc1dt = (c0 * Ca_cyt * kf) + (-c1*kb)
    J_cytB = (-Ca_cyt * c0 * kf) + (c1 * kb)

    # flux from SERCA
    V_max_serca = 250
    ka = 0.1 # converted to uM; it was 100 nM in the manuscript
    n = 2
    J_serca = V_max_serca * ((Ca_cyt ** n) / ((Ca_cyt ** n) + (ka ** n)))

    # flux due to ER calcium leak
    k_leak = 0.2
    J_ER_leak = k_leak*(Ca_ER-Ca_cyt)

    # flux due to ER calcium buffer
    kf_erB = 1
    kb_erB = 1200
    B0 = z[]
    B1 = z[]
    dB0dt = (-B0 * Ca_ER * kf_erB) + B1 * kb_erB
    dB1dt = (-B1 * kb_erB) + (B0 * Ca_ER * kf_erB)
    J_ERB = -Ca_cyt * B0 * kf_erB + B1 * kb_erB

#IP3 Dynamics
    # flux due to IP3 kinase
    V_max_ip3K = 5
    ka_ip3K1 = 0.4
    ka_ip3K2 = 10
    J_IP3_K = V_max_ip3K * (Ca_cyt**4)/((Ca_cyt**4)+(ka_ip3K1**4)) * IP3/(IP3 + ka_ip3K2)


