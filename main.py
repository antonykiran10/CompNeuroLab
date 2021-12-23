import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import time as clock

os.chdir('/home/antony/mycodes/anup_astron')
start_time = clock.time()

def model(z, t):
    # PMCA flux
    Ca_cyt = z[8]
    m0 = z[0]
    m1 = z[1]
    m2 = z[2]
    kf1 = 15
    kf2 = 20
    kf3 = 100
    k1 = 0.6
    kb1 = 20
    J_pmca = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m0 * k1)
    dm0dt = (-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m2 * kf3)
    dm1dt = (-m1 * kb1) + (m0 * Ca_cyt * kf1) + (-m1 * kf2)
    dm2dt = (-m2 * kf3) + (m1 * kf2)

    # cytosolic buffer
    c0 = z[3]
    c1 = z[4]
    kf = 60
    kb = 1200
    J_cytB = (-Ca_cyt * c0 * kf) + c1 * kb
    dc0dt = (-c0 * Ca_cyt * kf) + (c1 * kb)
    dc1dt = (-c0 * kb) + (c0 * Ca_cyt * kf)

    # SERCA FLUX
    V_max_serca = 250
    ka = 0.1
    n = 2
    J_serca = V_max_serca * ((Ca_cyt ** n) / ((Ca_cyt ** n) + (ka ** n)))

    # ER leak flux
    k_leak = 0.2
    Ca_ER = z[9]
    J_ERleak = k_leak * (Ca_ER - Ca_cyt)

    # ER buffer flux
    kfl = 1
    kbl = 80
    B0 = z[5]
    B1 = z[6]
    J_ERB = -Ca_cyt * B0 * kf + B1 * kbl
    dB0dt = (-B0 * Ca_ER * kfl) + B1 * kbl
    dB1dt = (-B1 * kbl) + (B0 * Ca_ER * kfl)

    # IP3 flux
    h = z[7]
    V_max_ip3R = 0.65
    a1 = 400
    a2 = 0.2
    a3 = 400
    a5 = 20
    b1 = 52
    b2 = 0.2
    b3 = 377.36
    b5 = 1.64
    d1 = b1 / a1
    d2 = b2 / a2
    d3 = b3 / a3
    d5 = b5 / a5
    f_noise = 0.2 / 100  # scaling factor for noise
    IP3_conc = z[10]
    q2 = d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    n = IP3_conc / (IP3_conc + d1)
    m = Ca_cyt / (Ca_cyt + d5)
    alpha_h = (a2 * d2) * ((IP3_conc + d1) / (IP3_conc + d3))
    beta_h = a2 * Ca_cyt
    varience = (alpha_h * (1 - h) + beta_h * h) / n
    std_dev = math.sqrt(abs(varience))
    J_IP3R = V_max_ip3R * (m ** 3) * (n ** 3) * (h ** 3) * (Ca_ER - Ca_cyt)
    dhdt = alpha_h * (1 - h) - beta_h * h + f_noise * abs(np.random.normal(0, std_dev))

    # IP3 kinase
    V_max_ip3K = 5
    ka1 = 0.4
    ka2 = 10
    J_IP3K = V_max_ip3K * ((Ca_cyt ** 4) / ((Ca_cyt ** 4) + (ka1 ** 4))) * (IP3_conc / (IP3_conc + ka2))

    # IP3 phosphatase
    V_max_ip3P = 1.25
    IP3_base = 0.160
    J_IP3P = V_max_ip3P * (IP3_conc - IP3_base)

    # IP3 PLC_delta
    V_max_ip3PLC = 0.2
    ka1_plc = 1
    ka2_plc = 1.5
    J_IP3PLC = (V_max_ip3PLC / (1 + (IP3_conc / ka1_plc))) * ((Ca_cyt ** 2) / ((Ca_cyt ** 2) + ka2_plc))

    # mGluR flux
    kGlu = 160
    Glu_conc = z[11]
    if 100 < t < 160:
        Glu_conc = 100
    # else:
    #     dGludt = (-Glu_conc * kGlu)
    dGludt = (-Glu_conc * kGlu)
    V_max_mGluR = 0.65
    ka_glu = 11
    J_mGluR = V_max_mGluR * ((Glu_conc ** 2) / ((Glu_conc ** 2) + (ka_glu ** 2)))

    volume_ratio = 2

    # final model

    dca_cytdt = J_pmca + J_cytB - J_serca + J_ERleak + J_IP3R
    dca_erdt = (J_serca - J_ERleak - J_IP3R) * volume_ratio + J_ERB
    dIP3dt = J_mGluR - J_IP3K - J_IP3P + J_IP3PLC

    dzdt = [dm0dt, dm1dt, dm2dt, dc0dt, dc1dt, dB0dt, dB1dt, dhdt, dca_cytdt, dca_erdt, dIP3dt, dGludt]
    return dzdt


initial_vals = [0, 0, 0, 0, 0, 0, 0, 0, 0.08, 400, 0.16, 0]  # [mo, m1, m2, c0, c1, B0, B1, h, ca_cyt, ca_er, IP3, Glu]

time_points = np.linspace(0, 200, int(200 / 0.00005))

result = odeint(model, initial_vals, time_points)

mydata = pd.DataFrame({
    'time': time_points,
    'Cal_cyto': result[:, 8],
    'Cal_ER': result[:, 9],
    'Glu': result[:, 10]
})

print("--- %s seconds ---" % (clock.time() - start_time))

mydata.to_csv('anup.csv')
plt.plot(time_points, result[:, 8], label='Ca_cyto')
# plt.plot(time_points,result[:,9], label="Ca_ER")
plt.grid()
plt.legend()
plt.show()
plt.plot(time_points, result[:, 9], label="Ca_ER")
plt.grid()
plt.legend()
plt.show()
plt.plot(time_points, result[:, 10], label="Glu")
plt.grid()
plt.legend()
plt.show()
