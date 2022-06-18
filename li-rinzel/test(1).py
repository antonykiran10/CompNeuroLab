import numpy as np
import matplotlib.pyplot as plt

v1 = np.random.normal(3 * 6, 0.2)
v2 = np.random.normal(3 * 0.11, 0.2)
v3 = np.random.normal(3 * 0.9, 0.2)

k3 = 0.1
c1 = 0.185

Dca = 20
Der = 0.2

d1 = 0.13
d2 = 1.049
d3 = 0.9434
d5 = 0.08234

a2 = 0.2

dx = 1
dt = 0.01


def dcal_cyto(cal_cyto, cal_ER, h, ip3_conc):
    return (-(c1 * v1 * m_inf(ip3_conc, cal_cyto) ** 3 * h ** 3 + v2) * (cal_cyto - cal_ER) - (v3 * cal_cyto ** 2) / (
            k3 ** 2 + cal_cyto ** 2)) * dt


def dcal_ER(cal_cyto, cal_ER, h, ip3_conc):
    return ((v1 * m_inf(ip3_conc, cal_cyto) ** 3 * h ** 3 + v2) * (cal_cyto - cal_ER) + (v3 * cal_cyto ** 2) / (
            c1 * (k3 ** 2 + cal_cyto ** 2))) * dt


def dh(h, ip3_conc, cal_cyto):
    return ((h_inf(ip3_conc, cal_cyto) - h) / tau_h(cal_cyto, ip3_conc)) * dt


def tau_h(cal_cyto, ip3_conc):
    return 1 / (a2 * (Q2(ip3_conc) + cal_cyto))


def Q2(ip3_conc):
    return d2 * ((ip3_conc + d1) / ip3_conc + d3)


def h_inf(ip3_conc, cal_cyto):
    return Q2(ip3_conc) / (Q2(ip3_conc) + cal_cyto)


def m_inf(ip3_conc, cal_cyto):
    return (ip3_conc / (ip3_conc + d1)) * (cal_cyto / (cal_cyto + d5))


N = int(100 / 0.001)
cal_cytoVal = np.arange(0, 100, 0.001)
print (cal_cytoVal)
m_infVal = np.zeros(N, float)
h_infVal = np.zeros(N, float)
# p_open = np.zeros(N, float)

for i in range(0, N):
    m_infVal[i] = m_inf(0.5, cal_cytoVal[i])
    h_infVal[i] = h_inf(0.5, cal_cytoVal[i])
    # p_open[i] = m_infVal[i] + h_infVal[i]


plt.plot(cal_cytoVal, h_infVal, label="h")
plt.plot(cal_cytoVal, m_infVal, label="m")
# plt.plot(cal_cytoVal, p_open, label="p_open")
plt.legend()
plt.show()
