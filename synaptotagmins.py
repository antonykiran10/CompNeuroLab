import numpy as np
import matplotlib.pyplot as plt

syt4 = np.zeros(100, int)
s0_count = 100
s1_count = 0
s2_count = 0
kf = 153
kb = 3500
Ca_cyt = 100
pf = (kf * Ca_cyt) / ((3 * kf * Ca_cyt) + (3 * kb))
pb = kb / ((3 * kf * Ca_cyt) + (3 * kb))

print(pf, pb)
print(3*pf + 3 * pb)

for j in range(0, 10000):

    for i in range(0, len(syt4)):
        if syt4[i] == 0 and np.random.binomial(1, 2 * pf):
            s0_count -= 1
            s1_count += 1
            syt4[i] = 1
        if syt4[i] == 1 and np.random.binomial(1, pf):
            s1_count -= 1
            s2_count += 1
            syt4[i] = 2
        if syt4[i] == 2 and np.random.binomial(1, 2 * pb):
            s2_count -= 1
            s1_count += 1
            syt4[i] = 1
        if syt4[i] == 1 and np.random.binomial(1, pb):
            s1_count -= 1
            s0_count += 1
            syt4[i] = 0

s0_count = np.count_nonzero(syt4 == 0)
s1_count = np.count_nonzero(syt4 == 1)
s2_count = np.count_nonzero(syt4 == 2)

print(s0_count, s1_count, s2_count)

print(syt4)
plt.plot(syt4)
plt.show()
