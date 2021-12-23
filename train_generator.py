import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

time_line = np.linspace(0, 200, 200)
wave = signal.square(2 * np.pi * 20 * time_line)
plt.plot(wave)
plt.show()


def input_train(t):
    if 100 < t < 150 and wave[t] > 0:
        return wave[t]
    else:
        return 0


a = np.zeros(300)
for i in range(0, 300):
    a[i] = input_train(i)

plt.plot(a)
plt.show()
