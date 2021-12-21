import numpy as np
import matplotlib.pyplot as plt

syt7 = np.zeros(100, int)

s7_0_count = 100
s7_1_count = 0
s7_2_count = 0
s7_3_count = 0
s7_4_count = 0
s7_5_count = 0

kf7 = 3.82
kb7 = 60
Ca_cyt = 100
pf7 = (kf7 * Ca_cyt) / ((6 * kf7 * Ca_cyt) + (6 * kb7))
pb7 = kb7 / ((6 * kf7 * Ca_cyt) + (6 * kb7))

print(pf7, pb7)
print(6 * pf7 + 6 * pb7)

for j in range(0, 10000):

    for i in range(0, len(syt7)):
        if syt7[i] == 0 and np.random.binomial(1, 5 * pf7):
            s7_0_count -= 1
            s7_1_count += 1
            syt7[i] = 1
        if syt7[i] == 1 and np.random.binomial(1, 4*pf7):
            s7_1_count -= 1
            s7_2_count += 1
            syt7[i] = 2
        if syt7[i] == 2 and np.random.binomial(1, 3 * pf7):
            s7_2_count -= 1
            s7_3_count += 1
            syt7[i] = 3
        if syt7[i] == 3 and np.random.binomial(1, 2*pf7):
            s7_3_count -= 1
            s7_4_count += 1
            syt7[i] = 4
        if syt7[i] == 4 and np.random.binomial(1, pf7):
            s7_4_count -= 1
            s7_5_count += 1
            syt7[i] = 5

        if syt7[i] == 5 and np.random.binomial(1, 5 * pb7):
            s7_5_count -= 1
            s7_4_count += 1
            syt7[i] = 4
        if syt7[i] == 4 and np.random.binomial(1, 4*pb7):
            s7_4_count -= 1
            s7_3_count += 1
            syt7[i] = 3
        if syt7[i] == 3 and np.random.binomial(1, 3 * pb7):
            s7_3_count -= 1
            s7_2_count += 1
            syt7[i] = 2
        if syt7[i] == 2 and np.random.binomial(1, 2*pb7):
            s7_2_count -= 1
            s7_1_count += 1
            syt7[i] = 1
        if syt7[i] == 1 and np.random.binomial(1, pb7):
            s7_1_count -= 1
            s7_0_count += 1
            syt7[i] = 0

s7_0_count = np.count_nonzero(syt7 == 0)
s7_1_count = np.count_nonzero(syt7 == 1)
s7_2_count = np.count_nonzero(syt7 == 2)
s7_3_count = np.count_nonzero(syt7 == 3)
s7_4_count = np.count_nonzero(syt7 == 4)
s7_5_count = np.count_nonzero(syt7 == 5)

print(s7_0_count, s7_1_count, s7_2_count, s7_3_count, s7_4_count, s7_5_count)

print(syt7)
plt.plot(syt7)
plt.show()
