import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def smpl_tcell_net(t,y):
    sol = np.zeros(y.shape)
    a = 1.; b = 1.; c = 1.; d = 1.; e = 1.; f = 1.; g = 1.; h = 1.; i = 1.; j = 1.

    sol[0] = a*y[3] +  b*y[0]*(0.5 - c*y[2])
    sol[1] = d*y[2] - e*y[1]
    sol[2] = f*y[3] - g*y[1]*y[2]
    sol[3] = h*y[4]*(i - y[3]) - y[3]
    sol[4] = -j*y[4]*y[0]

    return sol

def anton_tcell_net(t,y):
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

sol = solve_ivp(anton_tcell_net, [0.001,20.], np.asarray([0.,0.,0.,0.,1.]))
print(np.size(sol.t))
plt.plot(sol.t,sol.y[0], sol.t,sol.y[1], sol.t,sol.y[2], sol.t,sol.y[3], sol.t,sol.y[4])
#plt.ylim([0,100])
plt.yscale('log')
plt.legend(('X', 'Y', 'Z', 'ac', 'ag'))
plt.show()
