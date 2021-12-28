import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def model(z, t):
    Glu = z
    k_Glu = 160
    dGludt = (-Glu * k_Glu)
    return dGludt


initial_vals = [200]

time_points = np.linspace(0, 10, int(10 / 0.001))

result = odeint(model, initial_vals, time_points)

plt.plot(time_points,result, label="Glu")
plt.grid()
plt.legend()
plt.show()

