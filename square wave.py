from scipy import signal

import matplotlib.pyplot as plot

import numpy as np

# Sampling rate 1000 hz / second

t = np.linspace(0, 200, 200)

plot.plot(t, signal.square(2 * np.pi *10* t))
plot.grid()
plot.show()