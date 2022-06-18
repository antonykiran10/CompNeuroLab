from matplotlib import pyplot as pl
import numpy as np

# constant
c0 = 2.0
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
# time
dt = 0.0001
time = 100
N = int(time / dt)

cCa = np.zeros(N, float)
cCa[0] = 0.03
h = np.zeros(N, float)
h[0] = 0.7

for i in range(1, N):
    if 40 < dt * i < 60:
        cIP3 = 0.5
    else:
        cIP3 = 0
    m_infinity = cIP3 / (cIP3 + d1)
    n_infinity = (cCa[i - 1]) / (cCa[i - 1] + d5)
    alpha_h = a2 * d2 * (cIP3 + d1) / (cIP3 + d3)
    beta_h = a2 * cCa[i - 1]
    ERcCa = (c0 - cCa[i - 1]) / c1

    Jchannel = c1 * v1 * (m_infinity ** 3) * (n_infinity ** 3) * (h[i - 1] ** 3) * (cCa[i - 1] - ERcCa)
    Jpump = (v3 * cCa[i - 1] ** 2) / (k3 ** 2 + cCa[i - 1] ** 2)
    Jleak = c1 * v2 * (cCa[i - 1] - ERcCa)

    dCa = (-Jleak - Jchannel - Jpump) * dt
    dh = (alpha_h * (1 - h[i - 1]) - beta_h * (h[i - 1])) * dt

    cCa[i] = cCa[i - 1] + dCa
    h[i] = h[i - 1] + dh
    # print(cCa[i])

pl.plot(cCa)
pl.xlabel("time")
pl.ylabel("[cCa] micro molar")
pl.show()
