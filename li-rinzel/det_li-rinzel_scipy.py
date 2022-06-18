import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

# ip3_conc = []
def model(z, t):
    cal_cyto = z[0]
    h = z[1]
    c1 = 0.185
    v1 = 6.0
    v2 = 0.11
    v3 = 0.9
    k3 = 0.1
    d1 = 0.13
    d2 = 1.049
    d3 = 0.9434
    d5 = 0.08234
    a2 = 0.2
    c0 = 2.0
    if t > 150 and t < 160:
        IP3_conc = 1
    else:
        IP3_conc = 0

    # ip3_conc.append(IP3_conc)
    cal_ER = (c0 - cal_cyto) / c1
    m_inf = IP3_conc / (IP3_conc + d1)
    n_inf = cal_cyto / (cal_cyto + d5)
    alpha_h = a2 * d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    beta_h = a2 * cal_cyto
    J_channel = (c1 * v1 * (m_inf ** 3) * (n_inf ** 3) * (h**3) ) * (cal_cyto - cal_ER)
    J_pump = (v3 * (cal_cyto ** 2)) / ((k3 ** 2) + (cal_cyto ** 2))
    J_leak = c1 * v2 * (cal_cyto - cal_ER)
    dCadt = -J_channel - J_pump - J_leak
    dhdt = alpha_h*(1 - h) - beta_h * h
    dzdt = [dCadt, dhdt]
    return dzdt


initial_vals = [0.0, 0]

time_points = np.linspace(0, 300, int(100 / 0.001))
# ip3_conc = np.full(int(100/0.001), 0.5)
# ip3_conc = 0

result = odeint(model, initial_vals, time_points)
print(result[:,0])
# mydata= pd.DataFrame({
#     'time': time_points,
#     'Cal_cyto' = result[1]
# })
plt.plot(time_points,result[:,0])
plt.plot(time_points,result[:,1])
# plt.plot(time_points,ip3_conc)
plt.xlabel('Time (s)')
plt.ylabel('Calcium concentration (micro Molar)')
plt.title('Calcium concentration v/s Time')
plt.show()