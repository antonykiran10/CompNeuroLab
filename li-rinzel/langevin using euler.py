import numpy as np
import time as clock
import pandas as pd
import matplotlib.pyplot as plt

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

t = 200
dt = 0.01
N = int(t / dt)

start_time = clock.time()
#pool = np.random.normal(0, 0.1, 1000)


def dcal_cyto(cal_cyto, IP3_conc, h):
    cal_ER = (c0 - cal_cyto) / c1
    m_inf = IP3_conc / (IP3_conc + d1)
    n_inf = cal_cyto / (cal_cyto + d5)
    J_channel = (c1 * v1 * (m_inf ** 3) * (n_inf ** 3) * (h ** 3)) * (cal_cyto - cal_ER)
    J_pump = (v3 * (cal_cyto ** 2)) / ((k3 ** 2) + (cal_cyto ** 2))
    J_leak = c1 * v2 * (cal_cyto - cal_ER)
    # print(-J_channel , J_pump , J_leak)
    return (-J_channel - J_pump - J_leak) * dt


def dh(cal_cyto, IP3_conc, h):
    alpha_h = a2 * d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    beta_h = a2 * cal_cyto
    #print(alpha_h,cal_cyto)
    n = 5
    varience = (alpha_h * (1 - h) + beta_h * h) / n
    std_dev = np.sqrt(varience*dt*2)
    #print(varience)
    value_of_h = ((alpha_h * (1 - h) - beta_h * h))*dt + np.random.normal(loc=0.0, scale=std_dev)*dt

    # while True:
    #     print("inside")
    #     value_of_h = abs((alpha_h * (1 - h) - beta_h * h ) )*dt + np.random.standard_normal(1)
    #     if 0 <= value_of_h <= 1:
    #         break

    return value_of_h


def IP3_input(t,stim):
    if 100<t:
        return stim/100
    else:
        return 0


def my_looper(stimulus):
    stimulus = stimulus/100
    calcium_concentration = np.zeros(N, float)
    h_gate = np.zeros(N, float)
    time = np.zeros(N, float)
    calcium_concentration[0] = 0.035
    h_gate[0] = 0.7

    for i in range(1, N):
        calcium_concentration[i] = calcium_concentration[i - 1] + dcal_cyto(calcium_concentration[i - 1], IP3_input(i * dt,stimulus),
                                                                            h_gate[i - 1])
        #h_gate[i] = h_gate[i - 1] + dh(calcium_concentration[i - 1], IP3_input(i * dt), h_gate[i - 1])
        while True:
            h_gate[i] = h_gate[i - 1] + dh(calcium_concentration[i - 1], IP3_input(i * dt, stimulus), h_gate[i - 1])
            if 0 < h_gate[i] < 1:
                break
        time[i] = time[i - 1] + dt
    # df = pd.DataFrame({"Time (ms)": time, "Calcium": calcium_concentration, "h": h_gate})
    # df.to_csv("/home/antony/mycodes/li-rinzel/data_langevin/output_" + str(stimulus) + "_current" + ".csv", index=False)
    plt.plot(time, calcium_concentration, label="Ca in cyto")
    # # plt.plot(time, h_gate, label="H gate")
    # plt.grid()
    # plt.legend()
    # plt.show()
        #print(i)
for i in range(0,100,1):
    print(i)
    my_looper(i)

print("--- %s seconds ---" % (clock.time() - start_time))
# plt.plot(time, calcium_concentration, label="Ca in cyto")
# # # plt.plot(time, h_gate, label="H gate")
# plt.grid()
# plt.legend()
# plt.show()
