from ising_xy import IsingModelXY
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

ising_8 = IsingModelXY(1, 8, 1, 5000)
Ts = np.linspace(1,4,100)
ising_8.temperature_sweep(Ts)

with open('%s.csv' %'b_8', 'r') as my_file:
    m= np.loadtxt(my_file)

plt.figure()
plt.plot(Ts, m)
plt.show()
