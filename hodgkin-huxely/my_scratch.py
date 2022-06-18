import numpy as np
import matplotlib.pyplot as plt
import time as clock
import pandas as pd

start_time = clock.time()

# Vr=-65 # -65 mv resting membrane potential

T = 100
dt = 0.01
N = int(T / dt)


# input current
# unit μA/cm^2

def input_current(t):
    if 0 <= t < 5:
        return 3
    elif 11 <= t < 20:
        return 3
    else:
        return 0


# Value of the membrane Capacitance
# unit μF/cm^2

Cm = 1

# Value for reversal potentials
# unit : milli Volt (mV)

ENa = 115
EK = -12
EL = 10.6

# Value for conductance
# unit mS/cm^2

gNa = 120
gK = 36
gL = 0.3


def alpha_n(u):
    return (0.1 - 0.01 * u) / (np.exp(1 - 0.1 * u) - 1)


def beta_n(u):
    return 0.125 * np.exp(-u / 80)


def alpha_m(u):
    return (2.5 - 0.1 * u) / (np.exp(2.5 - 0.1 * u) - 1)


def beta_m(u):
    return 4 * np.exp(-u / 18)


def alpha_h(u):
    return 0.07 * np.exp(-u / 20)


def beta_h(u):
    return 1 / (np.exp(3 - 0.1 * u) + 1)


def dx(x, alpha, beta):
    return (alpha * (1 - x) - (beta * x)) * dt


def du(u, m_, n_, h_, t):
    return (input_current(t) / Cm - (
            (u - ENa) * gNa * m_ * m_ * m_ * h_ + (u - EK) * gK * n_ * n_ * n_ * n_ + (u - EL) * gL) / Cm) * dt


membrane_voltage = np.zeros(N, float)
m = np.zeros(N, float)
n = np.zeros(N, float)
h = np.zeros(N, float)
time = np.zeros(N, float)
input_ = np.zeros(N, float)

membrane_voltage[0] = 0
m[0] = (alpha_m(membrane_voltage[0]) / (alpha_m(membrane_voltage[0]) + beta_m(membrane_voltage[0])))
n[0] = (alpha_n(membrane_voltage[0]) / (alpha_n(membrane_voltage[0]) + beta_n(membrane_voltage[0])))
h[0] = (alpha_h(membrane_voltage[0]) / (alpha_h(membrane_voltage[0]) + beta_h(membrane_voltage[0])))
time[0] = 0
input_[0] = 0

for i in range(1, N):
    print("loop ", i)
    print("Potential ", membrane_voltage[i - 1])
    print(time[i - 1])
    time[i] = time[i - 1] + dt
    n[i] = n[i - 1] + dx(n[i - 1], alpha_n(membrane_voltage[i - 1]), beta_n(membrane_voltage[i - 1]))
    m[i] = m[i - 1] + dx(m[i - 1], alpha_m(membrane_voltage[i - 1]), beta_m(membrane_voltage[i - 1]))
    h[i] = h[i - 1] + dx(h[i - 1], alpha_h(membrane_voltage[i - 1]), beta_h(membrane_voltage[i - 1]))
    membrane_voltage[i] = membrane_voltage[i - 1] + du(membrane_voltage[i - 1], m[i - 1], n[i - 1], h[i - 1], time[i])
    input_[i] = input_current(time[i - 1])

print("--- %s seconds ---" % (clock.time() - start_time))

# df = pd.DataFrame({"Time (ms)": time, "Membrane Voltage (mV)": membrane_voltage, "m": m, "n": n, "h": h})
# df.to_csv(".\\run_data\\output_"+str(start_time)+".csv", index=False)

fig, ax = plt.subplots()
ax.plot(time, membrane_voltage)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Voltage (mV)')
ax.set_title('Membrane potential v/s Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, input_)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Current (mA)')
ax.set_title('Input current v/s Time')
plt.grid()

fig, ax = plt.subplots()
ax.plot(time, m, label='m')
ax.plot(time, n, label='n')
ax.plot(time, h, label='h')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Gating variable')
ax.set_title('Gating Variables w.r.t Time')
plt.legend()
plt.grid()

plt.show()
