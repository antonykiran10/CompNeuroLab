import numpy as np
import matplotlib.pyplot as plt

# assuming w = 1 

def a(pos):
    return -pos

dt =  0.01
T = 10
N = int(T/dt)

t = np.zeros(N, float)
x = np.zeros(N, float)
v = np.zeros(N, float)

t[0] = 0
x[0] = 10
v[0] = 0

for i in range (1,N):
    t[i] = t[i-1] + dt
    x[i] = x[i-1] + v[i-1]*dt
    v[i] = v[i-1] + a(x[i])*dt

plt . plot (t, x, label =" Position ")
plt . plot (t, v, label =" Velocity ")
plt . plot (t,-x, label =" Acceleration ")
plt . xlabel (" Time ")
plt . legend ()
plt . show ()