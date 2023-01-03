import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from ising_xy import IsingModelXY

T = np.arange(1.0, 4.05, 0.05)

with open('chis_8.csv') as f:
    b8 = np.loadtxt(f)

with open('chis_10.csv') as f:
    b10 = np.loadtxt(f)

with open('chis_12.csv') as f:
    b12 = np.loadtxt(f)

with open('chis_14.csv') as f:
    b14 = np.loadtxt(f)

with open('chis_16.csv') as f:
    b16 = np.loadtxt(f)

Tc = 2.18
L = [8, 10, 12, 14, 16]
nu = 0.68

for alpha in np.arange(-0.5,0.5,0.01):
    x_col = np.array([L[0]**(1/nu)*(T-Tc), L[1]**(1/nu)*(T-Tc), L[2]**(1/nu)*(T-Tc), L[3]**(1/nu)*(T-Tc), L[4]**(1/nu)*(T-Tc)]).flatten()
    idx = np.where((x_col < 10) & (x_col > -10))
    if len(x_col[idx]) > 0:
        y_col = np.array([b8*L[0]**(-alpha/nu), b10*L[1]**(-alpha/nu), b12*L[2]**(-alpha/nu), b14*L[3]**(-alpha/nu), b16*L[4]**(-alpha/nu)]).flatten()
        p = np.polyfit(x_col[idx], y_col[idx],2)
        chi_squared = np.sum((np.polyval(p, x_col[idx]) - y_col[idx]) ** 2)/np.std(y_col[idx])
        print(alpha, chi_squared)

alpha = -0.5

ax = plt.subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(direction='in')
plt.plot(L[0]**(1/nu)*(T-Tc), b8*L[0]**(-alpha/nu), '.-' , label='L=8')
plt.plot(L[1]**(1/nu)*(T-Tc), b10*L[1]**(-alpha/nu), '.-' , label='L=10')
plt.plot(L[2]**(1/nu)*(T-Tc), b12*L[2]**(-alpha/nu), '.-' , label='L=12')
plt.plot(L[3]**(1/nu)*(T-Tc), b14*L[3]**(-alpha/nu), '.-' , label='L=14')
plt.plot(L[4]**(1/nu)*(T-Tc), b16*L[4]**(-alpha/nu), '.-' , label='L=16')
plt.plot(np.repeat(0,100), np.linspace(0,100,100), 'k--')
plt.legend()
plt.xlim(-10,10)
plt.ylim(0,100)
plt.xlabel('$L^{1/\\nu}(T - T_c)$')
plt.ylabel('$mL^{-\\alpha/\\nu}$')
plt.savefig('alpha_-0.5.png', dpi=300)
plt.show()