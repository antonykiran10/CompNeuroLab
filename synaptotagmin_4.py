import numpy as np
import matplotlib.pyplot as plt

syt4 = np.zeros(100, int)
s4_0_count = 100
s4_1_count = 0
s4_2_count = 0
kf4 = 153
kb4 = 3500
Ca_cyt = 100
pf4 = (kf4 * Ca_cyt) / ((3 * kf4 * Ca_cyt) + (3 * kb4))
pb4 = kb4 / ((3 * kf4 * Ca_cyt) + (3 * kb4))

print(pf4, pb4)
print(3 * pf4 + 3 * pb4)

for j in range(0, 10000):

    for i in range(0, len(syt4)):
        if syt4[i] == 0 and np.random.binomial(1, 2 * pf4):
            s4_0_count -= 1
            s4_1_count += 1
            syt4[i] = 1
        if syt4[i] == 1 and np.random.binomial(1, pf4):
            s4_1_count -= 1
            s4_2_count += 1
            syt4[i] = 2
        if syt4[i] == 2 and np.random.binomial(1, 2 * pb4):
            s4_2_count -= 1
            s4_1_count += 1
            syt4[i] = 1
        if syt4[i] == 1 and np.random.binomial(1, pb4):
            s4_1_count -= 1
            s4_0_count += 1
            syt4[i] = 0

s4_0_count = np.count_nonzero(syt4 == 0)
s4_1_count = np.count_nonzero(syt4 == 1)
s4_2_count = np.count_nonzero(syt4 == 2)

print(s4_0_count, s4_1_count, s4_2_count)

print(syt4)
plt.plot(syt4)
plt.show()
