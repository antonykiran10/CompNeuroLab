import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
os.chdir('/home/antony/mycodes/li-rinzel/data')

run = time.strftime("%Y%m%d-%H%M%S")

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
    if t > 60 and t < 90:
        IP3_conc = 0.5
    else:
        IP3_conc = 0
    cal_ER = (c0 - cal_cyto) / c1
    m_inf = IP3_conc / (IP3_conc + d1)
    n_inf = cal_cyto / (cal_cyto + d5)
    alpha_h = a2 * d2 * ((IP3_conc + d1) / (IP3_conc + d3))
    beta_h = a2 * cal_cyto
    J_channel = (c1 * v1 * (m_inf ** 3) * (n_inf ** 3) * (h**3) ) * (cal_cyto - cal_ER)
    J_pump = (v3 * (cal_cyto ** 2)) / ((k3 ** 2) + (cal_cyto ** 2))
    J_leak = c1 * v2 * (cal_cyto - cal_ER)
    dCadt = -J_channel - J_pump - J_leak
    #dhdt = alpha_h*(1 - h) - beta_h * h
    a = (alpha_h * (1 - h)) - (beta_h * h)
    dhdt=a*(1+1/100*np.random.normal()) # 1 percent noise
    dzdt = [dCadt, dhdt]
    return dzdt

def loop(num):
    initial_vals = [0.035, 0.6]

    time_points = np.linspace(0, 100, int(100 / 0.001))

    result = odeint(model, initial_vals, time_points)
    # print(result[:,0])
    mydata= pd.DataFrame({
        'time': time_points,
        'Cal_cyto' : result[:,0],
        'h_gate' : result[:,1]
    })

    for i in range (0, len(mydata['h_gate'])):
        if mydata['h_gate'][i] > 1 or mydata['h_gate'][i] < 0:
            print ('h value error in run: ', i)
            continue
    print(result[:,0].max())

    mydata.to_csv(str(run)+'_run_'+str(num)+'.csv')
    # plt.plot(time_points,result[:,0])
    # plt.plot(time_points,result[:,1])
    # plt.show()
    return result[:,0].max()

calcium_peaks = []
run_id = []
for i in range(0,250):
    val = loop(i)
    calcium_peaks.append(val)
    run_id.append(i)

peak_data= pd.DataFrame({
    'run_id': run_id,
    'Cal_peak' : calcium_peaks
})

peak_data.to_csv(str(run)+'_summary.csv')
