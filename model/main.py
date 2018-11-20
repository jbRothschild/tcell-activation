import network as nk
import numpy as np
import matplotlib.pyplot as plt

def main():
    tmodel = nk.Models(np.asarray([0.,0.,0.,0.,100.]), nk.anton_net)
    tmodel.solve()
    print(np.size(tmodel.sol.t))
    plt.plot(tmodel.sol.t,tmodel.sol.y[0], tmodel.sol.t,tmodel.sol.y[1], tmodel.sol.t,tmodel.sol.y[2], tmodel.sol.t,tmodel.sol.y[3], tmodel.sol.t,tmodel.sol.y[4])
    #plt.ylim([0,100])
    plt.yscale('log')
    plt.legend(('X', 'Y', 'Z', 'ac', 'ag'))
    plt.show()

if __name__ == "__main__":
    main()
