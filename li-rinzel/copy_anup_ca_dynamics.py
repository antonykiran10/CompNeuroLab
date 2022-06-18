import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import math

os.chdir('/home/antony/mycodes/li-rinzel/data/anup')

run = time.strftime("%Y%m%d-%H%M%S")


def model(z, t):
    # def PMCA(z, t):
    m0 = z[0]
    m1 = z[1]
    m2 = z[2]
    Ca_cyt = z[8]
    Ca_ER = z[9]
    # IP3_conc = z[5]
    if t > 60 and t < 90:
        IP3_conc = 0.175
    else:
        IP3_conc = 0.16
    kf1 = 15
    kb1 = 20
    kf2 = 20
    kf3 = 100
    k1 = 0.6
    PMCA_density = 400

    dm0dt = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m2 * kf3)
    dm1dt = (-m1 * kb1) + (m0 * Ca_cyt * kf1) + (-m1 * kf2)
    dm2dt = (-m2 * kf3) + (m1 * kf2)
    J_pmca = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m0 * k1)

    # def cyt_buffer(z,t):
    c0 = z[5]
    c1 = z[6]
    kf = 60
    kb = 1200
    buffer_conc = 50

    dc0dt = (-c0 * Ca_cyt * kf) + (c1 * kb)
    dc1dt = (-c0 * kb) + (c0 * Ca_cyt * kf)
    J_cytB = (-Ca_cyt * c0 * kf) + c1 * kb

    # def SERCA(Ca_cyt):
    V_max_serca = 250
    ka = 0.1
    n = 2
    J_serca = V_max_serca * ((Ca_cyt ** 2) / ((Ca_cyt ** 2) + (ka ** 2)))

    # def ER_leak(Ca_ER, Ca_Cyto):
    k_leak = 0.2
    J_ERleak = k_leak * (Ca_ER - Ca_cyt)

    # IP3 receptor
    h = z[7]
    # IP3_conc = z[]
    a1 = 400
    a2 = 0.2
    a3 = 400
    a4 = 0.2
    a5 = 20
    b1 = 52
    b2 = 0.2
    b3 = 377.36
    b4 = 0.02
    b5 = 1.64
    d1 = b1 / a1
    d2 = b2 / a2
    d3 = b3 / a3
    # d4 = b4 / a4
    d5 = b1 / a5
    V_max_ip3R = 6.0
    n = 5
    n_inf = IP3_conc / (IP3_conc + d1)
    m_inf = Ca_cyt / (Ca_cyt + d5)
    alpha_h = a2 * d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    beta_h = a2 * Ca_cyt

    varience = (alpha_h * (1 - h) + beta_h * h) / n
    std_dev = math.sqrt(abs(varience))
    # print(std_dev)
    randomness = np.random.normal(0, std_dev)
    dhdt = alpha_h*(1 - h) - beta_h * h #+ (0.01/100)*randomness
    # print(dhdt)
    q2 = d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    J_ip3r = (V_max_ip3R * (m_inf ** 3) * (n_inf ** 3) * (h ** 3)) * (Ca_ER - Ca_cyt)

    # ER calcium buffer
    kfl = 1
    kbl = 80
    B0 = z[3]
    B1 = z[4]
    buffer_conc = 10

    dB0dt = (-B0 * Ca_ER * kfl) + B1 * kbl
    dB1dt = (-B1 * kbl) + (B0 * Ca_ER * kfl)
    J_ERB = -Ca_cyt * B0 * kf + B1 * kbl

    # # synaptogamins
    # kf4 = 153
    # kb4 = 3500
    # b4_syt = 0.25
    # kf7 = 3.82
    # kb7 = 60
    # J_syt4 = (-2 * s0 * kf4 * Ca_cyt) + (s1 * kb4) + (-s1 * kf4 * Ca_cyt) + (2 * s2 * kb4 * b4_syt)
    # J_syt7 = (-5 * Y0 * kf7 * Ca_cyt) + (Y1 * kb4) + (-4 * Y1 * kf7 * Ca_cyt) + (2 * Y2 * kb7 * b7) + (-3 * Y2 * kf7 * Ca_cyt) + (3 * Y3 * kf7 * Ca_cyt) + (-2 * Y3 * kf7 * Ca_cyt) + (4*Y4 * kb7 * (b7 ** 3)) + (-Y4 * kf7 * Ca_cyt) + (5 * Y5 * kb7 * (b7 ** 4))
    # J_syt = J_syt4 + J_syt7

    # final equations
    V_cyt = 0.37
    V_ER = 0.185
    dCal_cytdt = J_pmca + J_cytB - J_serca + J_ERleak + J_ip3r #+ J_syt
    dCa_ERdt = (J_serca - J_ERleak - J_ip3r)*(V_cyt/V_ER) + J_ERB

    # #IP3 fluxes
    # V_max_mGluR = 0.65
    # ka_mGluR = 6
    # J_mGluR = V_max_mGluR * ((Glu**2) / ((Glu**2)+(ka_mGluR**2)))

    # dzdt = [dm0dt, dm1dt, dm2dt, dCal_cytdt, dCa_ERdt, dc0dt, dc1dt, dhdt, dB0dt, dB1dt]
    dzdt = [dm0dt, dm1dt, dc0dt, dc1dt, dhdt, dB0dt, dB1dt, dm2dt, dCal_cytdt, dCa_ERdt]
    return dzdt

initial_vals = [0,0,0,0,0,0,0,0,0.08,400]

time_points = np.linspace(0, 100, int(100 / 0.001))

result = odeint(model, initial_vals, time_points)

mydata = pd.DataFrame({
    'time': time_points,
    'Cal_cyto': result[:, 8],
    'Cal_ER': result[:, 9]
})

# for i in range(0, len(mydata['h_gate'])):
#     if mydata['h_gate'][i] > 1 or mydata['h_gate'][i] < 0:
#         print('h value error in run: ', i)
#         continue
print(result[:, 4])

mydata.to_csv('anup.csv')
plt.plot(time_points,result[:,8], label='Ca_cyto')
plt.plot(time_points,result[:,9], label="Ca_ER")
plt.grid()
plt.legend()
plt.show()
