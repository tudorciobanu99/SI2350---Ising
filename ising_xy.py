import numpy as np
import matplotlib.pyplot as plt

class IsingModelXY:
    def __init__(self, L, J):
        self.L = L
        self.J = J
        self.N = L**3

    def reset_lattice(self):
        self.theta = np.random.uniform(0, 2*np.pi, (self.L,self.L,self.L))
    
    def energy_per_site(self, theta, i, j, k):
        e = -self.J*(np.cos(theta - self.theta[(i-1) % self.L, j, k])
        + np.cos(self.theta[(i+1) % self.L, j, k] - theta)
        + np.cos(theta - self.theta[i, (j-1) % self.L, k])
        + np.cos(self.theta[i, (j+1) % self.L, k] - theta)
        + np.cos(theta - self.theta[i, j, (k-1) % self.L])
        + np.cos(self.theta[i, j, (k+1) % self.L] - theta))
        return e

    def MCMC(self, T):
        for i in range(self.L):
            for j in range(self.L):
                for k in range(self.L):
                    e = self.energy_per_site(self.theta[i,j,k], i, j, k)
                    theta_candidate = np.random.uniform(0, 2*np.pi)
                    de = self.energy_per_site(theta_candidate,i,j,k) - e
                    if de <= 0 or (np.exp(-de/T) > np.random.uniform(0,1)):
                        self.theta[i,j,k] = theta_candidate

    def run_MCMC(self, T, num_warmup, num_sample):
        mabs = mabs2 = mabs4 = e = e2 = 0
        err = np.zeros((num_sample,5))
        for i in range(num_warmup):
            self.MCMC(T)
        for i in range(num_sample):
            self.MCMC(T)
            dm = np.abs(self.magnetization())
            de = self.compute_total_energy()
            mabs += dm
            e += de
            mabs2 += dm*dm
            mabs4 += dm*dm*dm*dm
            e2 += de*de
            err[i,:] = np.array([dm, dm*dm, dm*dm*dm*dm, e, e*e])
        mabs /= num_sample
        mabs2 /= num_sample
        mabs4 /= num_sample
        e /= num_sample
        e2 /= num_sample
        err = np.sqrt(1/(num_sample-1)*np.sum((err - np.array([mabs, mabs2, mabs4, e, e2]))**2, axis=0))/np.sqrt(num_sample)
        return mabs, mabs2, mabs4, e, e2, err

    def compute_total_energy(self):
        E = np.sum(-self.J*(np.cos(self.theta - np.roll(self.theta, 1, axis=0)) + np.cos(self.theta - np.roll(self.theta, -1, axis=0))
        + np.cos(self.theta - np.roll(self.theta, 1, axis=1)) + np.cos(self.theta - np.roll(self.theta, -1, axis=1))
        + np.cos(self.theta - np.roll(self.theta, 1, axis=2)) + np.cos(self.theta - np.roll(self.theta, -1, axis=2))))
        return E

    def magnetization(self):
        mx = np.sum(np.cos(self.theta))
        my = np.sum(np.sin(self.theta))
        return np.sqrt(mx**2 + my**2)/self.N

    def temperature_sweep(self, Ts, num_warmup, num_sample):
        mabs = np.zeros(len(Ts))
        chi = np.zeros(len(Ts))
        c = np.zeros(len(Ts))
        binder = np.zeros(len(Ts))
        err = np.zeros((len(Ts),5))
        for i in range(len(Ts)):
            T = Ts[i]
            self.reset_lattice()
            mabs[i], mabs2, mabs4, e, e2, err[i,:] = self.run_MCMC(T, num_warmup, num_sample)
            chi[i] = (mabs2 - mabs[i]*mabs[i])/T*self.N
            c[i] = (e2 - e*e)/(T*T*self.N)
            binder[i] = 1.0 - mabs4/(3*mabs2*mabs2)
            print(T, mabs[i], chi[i], c[i], binder[i])
        self.write_to_file('ms_%i'%self.L, mabs)
        self.write_to_file('chis_%i'%self.L, chi)
        self.write_to_file('cs_%i'%self.L, c)
        self.write_to_file('bs_%i'%self.L, binder)
        self.write_to_file('errs_%i'%self.L, err)
    
    def write_to_file(self, name, quantity):
        with open('%s.csv' %name, 'w') as my_file:
            np.savetxt(my_file, quantity)

    def equilibration_time(self, T, N):
        self.reset_lattice()
        e = np.zeros(N)
        for i in range(N):
            self.MCMC(T)
            e[i] = self.compute_total_energy()
        self.write_to_file('eq_%i'%self.L, e)




    
        
