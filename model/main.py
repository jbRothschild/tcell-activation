import network as nk
import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib.colors import LogNorm

def whole_plot_integration():
    """
    ############## WHOLE SYSTEM INTEGRATION ###########################
    for ag_0 in initial_dose:
        tmodel = nk.Models(np.asarray([0.,0.,0.,0.,ag_0]), nk.anton_net)
        tmodel.solve()

        fig = plt.figure()
        for i in np.arange(np.shape(tmodel.sol.y)[0]):
            plt.plot(tmodel.sol.t,tmodel.sol.y[i], color=plt.cm.Set2(i))
        plt.yscale('log'); plt.ylim([10**-2,10**5])
        plt.legend(tmodel.species)
        figname = "Network_ag_"+str(ag_0)
        plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    """
    return 0

def fold_expansion_antigen(initial_dose = np.logspace(0., 3., 10), DIR_OUTPUT='figures'):
    ## PLots similar to Figure 2 of Wilgreen et al. looking at fold expansion as a function of antigen concentration

    fig = plt.figure()
    fig.set_size_inches(10.2, 4.5)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for i, ag_0 in enumerate(initial_dose): # iterating over different initial antigen concentrations
        tmodel = nk.Models(np.asarray([0.,0.,0.,0.,ag_0]), nk.anton_net) # setting model
        tmodel.solve() #solving model
        time4 = np.argmax(tmodel.sol.t>1) # finding the right point in sol.t in which t closest to 2 time units (coulld change this)
        time7 = np.argmax(tmodel.sol.t>12)
        time9 = np.argmax(tmodel.sol.t>14)
        ax1.scatter(ag_0, tmodel.sol.y[0][time7]/tmodel.sol.y[0][time4], color=plt.cm.Set2(i)) # plotting the maximums as a function of the antigen concentration
        ax2.plot(tmodel.sol.t[time4:time9],tmodel.sol.y[0][time4:time9]/tmodel.sol.y[0][time4], color=plt.cm.Set2(i), label='ag='+str(int(ag_0))) # plotting
    ax1.set_yscale('log'); ax1.set_xscale('log'); ax1.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax1.set_xlabel(r'initial antigen concentration')
    ax2.set_yscale('log'); ax2.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax2.set_xlabel(r'time')
    ax2.legend();
    ax1.set_ylabel('Fold expansion (at time 12 normalized to concentration at time 1)')
    figname = "Fig2_Fold_expansion_vary_ag_0"
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    ## Similar plot
    fig = plt.figure()
    fig.set_size_inches(7.2, 4.5)
    ax1 = fig.add_subplot(111)
    for i, ag_0 in enumerate(initial_dose): # iterating over different initial antigen concentrations
        tmodel = nk.Models(np.asarray([0.,0.,0.,0.,ag_0]), nk.anton_net) # setting model
        tmodel.solve() #solving model
        ax1.plot(tmodel.sol.t,tmodel.sol.y[0], color=plt.cm.Set2(i),label='ag='+str(int(ag_0))) # plotting
    ax1.set_yscale('log'); ax1.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax1.set_xlabel(r'Time')
    ax1.legend();
    plt.ylabel('Concentration T')
    figname = "Network_vary_ag_0"
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    return 0

def fold_expansion_initialT(initial_T_dose = np.logspace(0., 6., 10), ag_0=200, DIR_OUTPUT='figures'):
    ## PLots similar to Figure 1 of Wilgreen et al. looking at fold expansion as a function of antigen concentration

    fig = plt.figure()
    fig.set_size_inches(10.2, 4.5)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    for i, tcell_0 in enumerate(initial_T_dose): # iterating over different initial antigen concentrations
        tmodel = nk.Models(np.asarray([tcell_0,0.,0.,0.,ag_0]), nk.anton_net) # setting model
        tmodel.solve() #solving model
        ax1.scatter(tcell_0, tmodel.sol.y[0][np.argmax(tmodel.sol.y[0])]/tcell_0, color=plt.cm.Set2(i), label=r'$T_{pr}(0)=$'+str(int(tcell_0))) # plotting the maximums as a function of the antigen concentration
        ax2.plot(tmodel.sol.t,tmodel.sol.y[0], color=plt.cm.Set2(i)) # plotting
    ax1.set_yscale('log'); ax1.set_xscale('log'); ax1.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax1.set_xlabel(r'initial concentration T cells')
    ax2.set_yscale('log');  ax2.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax2.set_xlabel(r'time')
    ax1.legend();
    ax1.set_ylabel('Fold expansion vary initial T cell'); ax2.set_ylabel("Concentration of proliferating T cells")
    plt.title('Antigen concentration at ag='+str(int(ag_0)))
    figname = "Fig1_Fold_expansion_vary_tcell_0_"+str(int(ag_0))
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')

    ## Similar plot
    fig = plt.figure()
    fig.set_size_inches(7.2, 4.5)
    ax1 = fig.add_subplot(111)
    for i, tcell_0 in enumerate(initial_T_dose): # iterating over different initial antigen concentrations
        tmodel = nk.Models(np.asarray([tcell_0,0.,0.,0.,ag_0]), nk.anton_net) # setting model
        tmodel.solve() #solving model
        ax1.plot(tmodel.sol.t,tmodel.sol.y[0]/tcell_0, color=plt.cm.Set2(i), label=r'$T_{pr}(0)=$'+str(int(tcell_0))) # plotting
    ax1.set_yscale('log'); ax1.set_ylim([10**-2,np.max(tmodel.sol.y[0])]); ax1.set_xlabel(r'Time')
    plt.ylabel('Concentration T'); plt.legend()
    plt.title('Antigen concentration at ag='+str(int(ag_0)))
    figname = "Network_vary_tcell_0_"+str(int(ag_0))
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf'); plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    return 0

def main():
    DIR_OUTPUT = "figures"
    if not os.path.exists(DIR_OUTPUT):
        os.makedirs(DIR_OUTPUT)
    tmodel = nk.Models(np.asarray([0.,0.,0.,0.,0.]), nk.anton_net)

    n = 10.; m = 20.
    initial_dose = np.logspace(0., 3., n); initial_T_dose = np.logspace(-2., 2, m)

    fold_expansion_antigen(initial_dose)
    fold_expansion_initialT(ag_0=30)
    fold_expansion_initialT(ag_0=200)
    fold_expansion_initialT(ag_0=1000)
    """
    ########### VARYING
    SOL_AG = []
    SOL_T = []

    for i in initial_dose:
        tmodel.change_init_conc(np.asarray([0.,0.,0.,0.,i]))
        tmodel.solve()
        SOL_AG.append(tmodel.sol.y[0])
        print(tmodel.sol.t[np.argmax(tmodel.sol.y[0])])

    SOL_array = (np.asarray(SOL_AG))
    SOL_array[SOL_array <= 0] = 10**(-7)
    maximum = int(np.log10(SOL_array).max())+1
    minimum = 0

    fig  = plt.figure()
    ax = fig.add_subplot(111)

    cax = ax.contourf(tmodel.sol.t, initial_dose, SOL_array, cmap=plt.cm.inferno, norm=LogNorm() , levels=np.logspace(minimum, maximum, 100))
    cbar = fig.colorbar(cax, ticks=[10**minimum, 10**int((maximum-minimum)/3), 10**int((maximum-minimum)*2/3), 10**maximum])
    for c in ax.collections:
        c.set_edgecolor("face")
    cbar.ax.set_ylabel(r'$Concentration\ of\ T_{cell}$')
    ax.set_title(r"")
    ax.set_xlabel(r"$Time$")
    ax.set_yscale('log')
    ax.set_ylabel(r"$ initial\ antigen\ concentration$")
    #ax.set_yscale("log")
    ax.minorticks_off()
    plt.savefig("T_cell_by_initial_ag.png")
    plt.close(fig)

    #################################


    for i in initial_dose:
        fig  = plt.figure()
        ax = fig.add_subplot(111)
        for j in initial_T_dose:
            print(j)
            tmodel.change_init_conc(np.asarray([j,0.,0.,0,i]))
            tmodel.solve()
            SOL_T.append(tmodel.sol.y[0])
            print(j)
        plt.plot(initial_T_dose, np.asarray(SOL_T).max(axis=1)/initial_T_dose)
        ax.set_title(r"")
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel(r"$Fold\ Expansion$")
        ax.set_xlabel(r"$Initial\ reg\ T\ cell\ concentration$")
        ax.minorticks_off()
        plt.savefig("T_cell_max_by_initial_"+str(i)+".png")
        plt.show()
        plt.close(fig)
        SOL_T = []

    SOL_T_array = (np.asarray(SOL_T))
    SOL_T_array[SOL_T_array <= 0] = 10**(-7)
    maximumT = int(np.log10(SOL_array).max())+1
    minimumT = 0


    fig  = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(initial_dose, [tmodel.sol.t[time] for time in np.argmax(SOL_array,axis=1)])
    ax.set_title(r"")
    ax.set_xscale('log')
    ax.set_ylabel(r"$Time of Maximum\ T\ cell\ concentration$")
    ax.set_xlabel(r"$ initial\ antigen\ concentration$")
    #ax.set_yscale("log")
    ax.minorticks_off()
    plt.savefig("T_cell_max_time_by_initial_ag.png")
    plt.show()
    """

if __name__ == "__main__":
    main()
