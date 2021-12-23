from scipy import signal

import matplotlib.pyplot as plot

import numpy as np

# Sampling rate 1000 hz / second

t = np.linspace(0, 200, 200)
a = signal.square(2 * np.pi *100* t)
print(np.count_nonzero(a == 1))
plot.plot(t, a)
plot.grid()
plot.show()