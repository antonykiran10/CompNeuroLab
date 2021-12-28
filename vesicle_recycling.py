import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib.pyplot as plt

# Recycling scheme
def recyle(z, t):
    docked = z[0]
    mobile = z[1]
    released = z[2]
    endocytosed = z[3]
    reacidified = z[4]
    # if 100 < t < 150:
    #     docked = 1
    r04 = 0.000417
    rCa4 = 115
    r07 = 0.000417
    rCa7 = 8
    kef = 0.66
    kaf = 0.6
    kes = 0.16
    kas = 0.052
    k_mob = 0.615
    k_doc = 0.75
    ddockeddt = k_doc * mobile - (r04 + rCa4) * docked
    dmobiledt = (k_mob * reacidified) - (k_doc * mobile) - (r07 + rCa7) * mobile
    dreleaseddt = (r04 + rCa4) * docked + (r07 + rCa7) * mobile - (kef + kes) * endocytosed
    dendodt = (kef + kes) * released - (kaf + kas) * endocytosed
    dreaciddt = (kaf + kas) * endocytosed - k_mob * reacidified
    # ddockeddt = k_doc * docked - (r04 + rCa4) * docked
    # dmobiledt = (k_mob * mobile) - (k_doc * mobile) - (r07 + rCa7) * mobile
    # dreleaseddt = (r04 + rCa4) * released + (r07 + rCa7) * released - (kef + kes) * released
    # dendodt = (kef + kes) * endocytosed - (kaf + kas) * endocytosed
    # dreaciddt = (kaf + kas) * reacidified - k_mob * reacidified
    dzdt = [ddockeddt, dmobiledt, dreleaseddt, dendodt, dreaciddt]
    return dzdt


docked_vesicles = 80
mobile_vesicles = 20
released_vesicles = 0
endocytosed_vesicles = 0
reacidified_vesicles = 0
initial_vals = [docked_vesicles, mobile_vesicles, released_vesicles, endocytosed_vesicles, reacidified_vesicles]

time_points = np.linspace(0, 200, int(200 / 0.00005))

result = odeint(recyle, initial_vals, time_points)
plt.plot(time_points, result)
print(result[:,0])
plt.show()
# print(result)

#
# mydata = pd.DataFrame({
#     'time': time_points,
#     'docked': result[:, 0],
#     'mobile': result[:, 1],
#     'released': result[:, 2],
#     'endocytosed': result[:, 3],
#     'reacidified': result[:, 4]
# })
#
# mydata.to_csv('vesicle_recycling.csv')
# plt.plot(time_points, result[:, 0], label='docked')
# plt.show()
