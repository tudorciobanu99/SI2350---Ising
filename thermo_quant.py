import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from ising_xy import IsingModelXY

T = np.arange(1.0, 4.05, 0.05)

with open('bs_8.csv') as f:
    b8 = np.loadtxt(f)

with open('bs_10.csv') as f:
    b10 = np.loadtxt(f)

with open('bs_12.csv') as f:
    b12 = np.loadtxt(f)

with open('bs_14.csv') as f:
    b14 = np.loadtxt(f)

with open('bs_16.csv') as f:
    b16 = np.loadtxt(f)

with open('ms_8.csv') as f:
    ms8 = np.loadtxt(f)

with open('ms_10.csv') as f:
    ms10 = np.loadtxt(f)

with open('ms_12.csv') as f:
    ms12 = np.loadtxt(f)

with open('ms_14.csv') as f:
    ms14 = np.loadtxt(f)

with open('ms_16.csv') as f:
    ms16 = np.loadtxt(f)

with open('errs_8.csv') as f:
    err8 = np.loadtxt(f)

with open('errs_10.csv') as f:
    err10 = np.loadtxt(f)

with open('errs_12.csv') as f:
    err12 = np.loadtxt(f)

with open('errs_14.csv') as f:  
    err14 = np.loadtxt(f)

with open('errs_16.csv') as f:
    err16 = np.loadtxt(f)

ax = plt.subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction='in')
plt.plot(T, b8,'.-', label='L=8')
plt.plot(T, b10,'.-', label='L=10')
plt.plot(T, b12,'.-', label='L=12')
plt.plot(T, b14,'.-', label='L=14')
plt.plot(T, b16,'.-', label='L=16')
plt.plot(T, np.repeat(1/3, len(T)), 'k--')
plt.plot(T, np.repeat(2/3, len(T)), 'k--')
plt.plot(np.repeat(2.18,len(T)), np.linspace(0,1,len(T)), 'k-')
plt.ylim(0.2, 0.8)
plt.xlim(1,4)
plt.xlabel('$T$')
plt.ylabel('$B$')
plt.legend()
plt.tight_layout()
plt.savefig('bs_zoom.png', dpi =  300)
plt.show()