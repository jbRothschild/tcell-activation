import numpy as np
from scipy.integrate import solve_ivp

def anton_net(t,y):
    """
    Args :
        t : time points
        y[0] : X ; y[1] : Y ; y[2] : Z ; y[3] : ac ; y[4] : ag ;
    """
    sol = np.zeros(y.shape)
    a = 1.3; b = 0.1; c = 0.1; d = 10.; e = 0.2; f = 100.; g = 0.1;
    sol[0] = y[0]*(1 - a/(1 + y[2])) + b*y[3]
    sol[1] = y[2]/(1 + y[2]) - c*y[1]
    sol[2] = y[3] - d*y[1]*y[2]/(1 + y[2])
    sol[3] = e*y[4]*(f - y[3]) - y[3]
    sol[4] = -g*y[4]*y[0]

    return sol

class Models:

    def __init__(self, initial_concentration, ntwk):
        self.init_conc = initial_concentration
        self.init_time = 0.001
        self.final_time = 20.
        self.network = ntwk
        self.time_eval = list(np.arange(self.init_time, self.final_time, (self.final_time-self.init_time)/1000.))
        self.species = (r'$T_{pr}$', r'$C_{IL2}$', r'$T_{reg}$', r'$T_{ac}$', r'$ag$')

    def change_time(self, begin, end):
        self.init_time = begin
        self.final_time = end

    def change_init_conc(self, new_conc):
        self.init_conc = new_conc

    def change_time_eval(self, intervals):
        self.time_eval = list(np.arange(self.init_time, self.final_time, (self.final_time-self.init_time)/intervals))

    def solve(self):
        self.sol = solve_ivp(self.network, [self.init_time, self.final_time], self.init_conc, t_eval=self.time_eval)
