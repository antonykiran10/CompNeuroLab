# Roth, B. J., et al. “A Mathematical Model of Agonist-Induced Propagation of Calcium Waves in Astrocytes.”
# Cell Calcium, vol. 17, no. 1, Jan. 1995, pp. 53–64. DOI.org (Crossref),
# https://doi.org/10.1016/0143-4160(95)90102-7.


import numpy as np
import time as clock
import matplotlib.pyplot as plt

start_time = clock.time()
v1 = 4 * 6
v2 = 4 * 0.011
v3 = 4 * 0.9

k3 = 0.1
c1 = 0.185

Dca = 20
Der = 0.2

d1 = 0.13
d2 = 1.049
d3 = 0.9434
d5 = 0.08234

a2 = 0.2

L = 80
T = 100
dx = 1
dt = 0.01
N = int(L / dx)
M = int(T / dt)


def ip3_input(t):
    # return 0.3
    if 60 <= t < 90:
        return 0.5
    else:
        return 0.02


def dcal_cyto(cal_cyto, cal_ER, h, ip3_conc):
    # return (-(c1 * (v1 * (m_inf(ip3_conc, cal_cyto) ** 3) * (h ** 3) + v2)) * (cal_cyto - cal_ER) -
    #         ((v3 * (cal_cyto ** 2)) / ((k3 ** 2) + (cal_cyto ** 2)))) * dt
    m = m_inf(ip3_conc, cal_cyto)
    return (-(c1 * ((v1 * m * m * m * h * h * h) + v2) * (cal_cyto - cal_ER)) - (
            (v3 * cal_cyto * cal_cyto) / (k3 * k3 + cal_cyto * cal_cyto))) * dt


def dcal_ER(cal_cyto, cal_ER, h, ip3_conc):
    # return ((v1 * (m_inf(ip3_conc, cal_cyto) ** 3) * (h ** 3) + v2) * (cal_cyto - cal_ER) + (
    #             (1 / c1) * ((v3 * (cal_cyto ** 2)) / (
    #             c1 * ((k3 ** 2) + (cal_cyto ** 2)))))) * dt
    m = m_inf(ip3_conc, cal_cyto)
    return ((((v1 * m * m * m * h * h * h) + v2) * (cal_cyto - cal_ER)) + ((1 / c1) * (
            (v3 * cal_cyto * cal_cyto) / (k3 * k3 + cal_cyto * cal_cyto)))) * dt


# def cal_diffusion_cyto(cal_cyto):


def dh(h, ip3_conc, cal_cyto):
    # return ((h_inf(ip3_conc, cal_cyto) - h) / tau_h(cal_cyto, ip3_conc)) * dt
    hinf = h_inf(ip3_conc, cal_cyto)
    tau = tau_h(cal_cyto, ip3_conc)
    return ((hinf - h) / tau) * dt


def tau_h(cal_cyto, ip3_conc):
    return 1 / (a2 * (Q2(ip3_conc) + cal_cyto))


def Q2(ip3_conc):
    # print (ip3_conc, d1, d3)
    # print(d2 * ((ip3_conc + d1) / ip3_conc + d3))
    return d2 * ((ip3_conc + d1) / (ip3_conc + d3))


def h_inf(ip3_conc, cal_cyto):
    q2 = Q2(ip3_conc)
    return q2 / (q2 + cal_cyto)


def m_inf(ip3_conc, cal_cyto):
    # print(ip3_conc, cal_cyto)
    # print(float((ip3_conc / (ip3_conc + d1)) * (cal_cyto / (cal_cyto + d5))))
    return (ip3_conc / (ip3_conc + d1)) * (cal_cyto / (cal_cyto + d5))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    calcium_concentration = np.zeros(M, float)
    calcium_concentration_ER = np.zeros(M, float)
    h_vals = np.zeros(M, float)
    time = np.zeros(M, float)
    ip3_concentration = np.zeros(M, float)
    m_infVals = np.zeros(M, float)
    h_infVals = np.zeros(M, float)

    calcium_concentration[0] = 0.035
    calcium_concentration_ER[0] = 50
    h_vals[0] = 0
    ip3_concentration[0] = 0.5
    m_infVals[0] = m_inf(ip3_concentration[0], calcium_concentration[0])
    h_infVals[0] = h_inf(ip3_concentration[0], calcium_concentration[0])

    for i in range(1, M):
        # print(calcium_concentration[i - 1])
        calcium_concentration[i] = calcium_concentration[i - 1] + dcal_cyto(calcium_concentration[i - 1],
                                                                            calcium_concentration_ER[i - 1],
                                                                            h_vals[i - 1], ip3_concentration[i - 1])
        calcium_concentration_ER[i] = calcium_concentration_ER[i - 1] + dcal_ER(calcium_concentration[i - 1],
                                                                                calcium_concentration_ER[i - 1],
                                                                                h_vals[i - 1], ip3_concentration[i - 1])
        h_vals[i] = h_vals[i - 1] + dh(h_vals[i - 1], ip3_concentration[i - 1], calcium_concentration[i - 1])
        print(h_vals[i], calcium_concentration[i], calcium_concentration_ER[i])
        time[i] = time[i - 1] + dt
        ip3_concentration[i] = ip3_input(time[i - 1])
        m_infVals[i] = m_inf(ip3_concentration[i], calcium_concentration[i])
        h_infVals[i] = h_inf(ip3_concentration[i], calcium_concentration[i])

    print("--- %s seconds ---" % (clock.time() - start_time))
    plt.plot(time, calcium_concentration, label="Ca in cyto")
    # plt.plot(time, calcium_concentration_ER, label="Ca in ER")
    # plt.plot(time, h_vals, label="h")
    # plt.plot(time, m_infVals, label="M_inf")
    # plt.plot(time, h_vals, label="H")
    plt.grid()
    plt.legend()
    plt.show()
