import numpy as np
import matplotlib.pyplot as plt


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


voltage = np.arange(-100, 200, 0.01)
print(len(voltage))
n = np.zeros(len(voltage), float)
m = np.zeros(len(voltage), float)
h = np.zeros(len(voltage), float)
n[0] = (alpha_n(voltage[0]) / (alpha_n(voltage[0]) + beta_n(voltage[0])))
m[0] = (alpha_m(voltage[0]) / (alpha_m(voltage[0]) + beta_m(voltage[0])))
h[0] = (alpha_h(voltage[0]) / (alpha_h(voltage[0]) + beta_h(voltage[0])))
dt = 0.0001
for i in range(1, len(voltage)):
    n[i] = n[i-1] + dx(n[i - 1], alpha_n(voltage[i - 1]), beta_n(voltage[i - 1]))
    m[i] = m[i-1] + dx(m[i - 1], alpha_m(voltage[i - 1]), beta_m(voltage[i - 1]))
    h[i] = h[i-1] + dx(h[i - 1], alpha_h(voltage[i - 1]), beta_h(voltage[i - 1]))

plt.plot(voltage, n)
plt.plot(voltage, m)
plt.plot(voltage, h)
plt.show()
