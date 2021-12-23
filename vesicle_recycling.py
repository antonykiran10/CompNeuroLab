import numpy as np
vesicles = np.zeros(1000)

Docked = 0
Mobile = 1
Released = 2
Endocytosed = 3
Reacidified = 4


# creaion of docked and mobile vesicle pool
fraction_of_mobile = 0.2
for i in range(0, len(vesicles)):
    if vesicles[i] == Docked and np.random.binomial(1, fraction_of_mobile):
        vesicles[i] = Mobile


Docked_count = np.count_nonzero(vesicles == Docked)
Mobile_count = np.count_nonzero(vesicles == Mobile)
Released_count = np.count_nonzero(vesicles == Released)
Reacidified_count = np.count_nonzero(vesicles == Reacidified)
# Recycling scheme
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
for i in range(0, len(vesicles)):
    if vesicles[i] == Docked and np.random.binomial(1,)

Docked_count = np.count_nonzero(vesicles == Docked)
Mobile_count = np.count_nonzero(vesicles == Mobile)

print( Docked_count, Mobile_count, Released_count, Reacidified_count)