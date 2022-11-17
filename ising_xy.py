import numpy as np
import matplotlib.pyplot as plt

class IsingModelXY:

    def __init__(self, T, L, J, num_steps):
        self.L = L
        self.T = T
        self.J = J
        self.num_steps = num_steps
        self.N = L**3
        self.theta = np.zeros((L, L, L))
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    self.theta[i, j, k] = np.random.uniform(0,2*np.pi)

    def find_neighbors(self, i, j, k):
        index = [self.L-1 if i == 0 else i-1, 0 if i == self.L-1 else i+1, self.L-1 if j == 0 else j-1, 0 if j == self.L-1 else j+1, self.L-1 if k == 0 else k-1, 0 if k == self.L-1 else k+1]
        return index
    
    def energy_per_site(self, theta, i, j, k):
        index = self.find_neighbors(i,j,k)
        E = -self.J*(np.cos(theta - self.theta[index[0],j,k])
        + np.cos(self.theta[index[1],j,k] - theta)
        + np.cos(theta- self.theta[i,index[2],k])
        + np.cos(self.theta[i,index[3],k] - theta)
        + np.cos(theta - self.theta[i,j,index[4]])
        + np.cos(self.theta[i,j,index[5]] - theta))
        return E
    
    def flip(self):
        i = np.random.randint(0,self.L)
        j = np.random.randint(0,self.L)
        k = np.random.randint(0,self.L)
        E = self.energy_per_site(self.theta[i,j,k], i, j, k)
        theta_candidate = np.random.uniform(0,2*np.pi)
        dE = self.energy_per_site(theta_candidate,i,j,k) - E
        if dE < 0 or (np.exp(-dE/self.T) > np.random.uniform(0,1)):
            self.theta[i,j,k] = theta_candidate

    def MCMC(self):
        for i in range(self.N):
            self.flip()

    def run_MCMC(self):
        m = 0
        m_sq = 0
        m_f = 0
        En = 0
        E_sq = 0
        for i in range(self.num_steps):
            if i > self.num_steps/10:
                mag = self.magnetization()
                m += np.abs(mag)
                m_sq += mag**2
                m_f += mag**4
                energy = self.compute_total_energy()
                En += energy
                E_sq += energy**2
            self.MCMC()
        m = m/self.num_steps
        m_sq = m_sq/self.num_steps
        m_f = m_f/self.num_steps
        En = En/self.num_steps
        E_sq = E_sq/self.num_steps
        return m, m_sq, m_f, En, E_sq

    def compute_total_energy(self):
        En = 0
        for i in range(self.L):
            for j in range(self.L):
                for k in range(self.L):
                    En += self.energy_per_site(self.theta[i, j, k], i, j, k)
        return En

    def magnetization(self):
        mx = np.sum(np.cos(self.theta))
        my = np.sum(np.sin(self.theta))
        return np.sqrt(mx**2 + my**2)/self.N

    def temperature_sweep(self, Ts):
        m_avg = np.zeros(len(Ts))
        susc = np.zeros(len(Ts))
        m_sq = np.zeros(len(Ts))
        m_f = np.zeros(len(Ts))
        En = np.zeros(len(Ts))
        E_sq = np.zeros(len(Ts))
        c = np.zeros(len(Ts))
        b = np.zeros(len(Ts))
        for i in range(len(Ts)):
            self.T = Ts[i]
            m_avg[i], m_sq[i], m_f[i], En[i], E_sq[i] = self.run_MCMC()
            susc[i] = 1/self.T*self.N*(m_sq[i] - m_avg[i]**2)
            c[i] = 1/(self.T**2*self.N)*(E_sq[i] - En[i]**2)
            b[i] = 1 - m_f[i]/(3*m_sq[i]**2)
        self.write_to_file('m_%i'%self.L, m_avg)
        self.write_to_file('chi_%i'%self.L, susc)
        self.write_to_file('c_%i'%self.L, c)
        self.write_to_file('b_%i'%self.L, b)
    
    def write_to_file(self, name, quantity):
        with open('%s.csv' %name, 'w') as my_file:
            np.savetxt(my_file, quantity)

    
        
