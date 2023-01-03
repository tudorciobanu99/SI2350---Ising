from ising_xy import IsingModelXY
import numpy as np
import matplotlib.pyplot as plt

with open('eq1_16.csv', 'r') as f:
    e1 = np.loadtxt(f)
with open('eq2_16.csv', 'r') as f:
    e2 = np.loadtxt(f)
with open('eq3_16.csv', 'r') as f:
    e3 = np.loadtxt(f)
with open('eq4_16.csv', 'r') as f:
    e4 = np.loadtxt(f)
with open('eq5_16.csv', 'r') as f:
    e5 = np.loadtxt(f)
with open('eq6_16.csv', 'r') as f:
    e6 = np.loadtxt(f)
with open('eq7_16.csv', 'r') as f:
    e7 = np.loadtxt(f)
with open('eq8_16.csv', 'r') as f:
    e8 = np.loadtxt(f)
with open('eq9_16.csv', 'r') as f:
    e9 = np.loadtxt(f)
plt.figure()
n = range(20000)
plt.plot(n, e1, color = 'red')
plt.plot(n, e2, color = 'yellow')
plt.plot(n, e3, color = 'orange')
plt.plot(n, e4, color = 'green')
plt.plot(n, e5, color = 'blue')
plt.plot(n, e6, color = 'purple')
plt.plot(n, e7, color = 'pink')
plt.plot(n, e8, color = 'brown')
plt.plot(n, e9, color = 'black')
plt.text(1000, -21000, "T = 1.0", ha="center", va="center", rotation=0, size=10,
    bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1))
plt.text(1000, -4500, "T = 4.0", ha="center", va="center", rotation=0, size=10,
    bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1))
plt.text(1000, -10000, "T = 2.2", ha="center", va="center", rotation=0, size=10,
    bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1))
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Energy')
plt.xlim(0,2000)
plt.tight_layout()
plt.ylim(-22500, -2500)
plt.savefig('eq_16.png', dpi =  300)
plt.show()