import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys

# Plot params
mpl.rcParams['figure.figsize'] = (10,7)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.linewidth'] = 2.5
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['font.serif'] = ['Times']

def make_axis(axe, logX, logY):
    if not logX:
        axe.xaxis.set_minor_locator(AutoMinorLocator())
    if not logY:
        axe.yaxis.set_minor_locator(AutoMinorLocator())
    axe.xaxis.set_tick_params(which='major', direction='in', size=10, width=1, top='on', bottom='on')
    axe.xaxis.set_tick_params(which='minor', direction='in', size=5, width=1, top='on', bottom='on')
    axe.yaxis.set_tick_params(which='major', direction='in', size=10, width=1, right='on', left='on')
    axe.yaxis.set_tick_params(which='minor', direction='in', size=5, width=1, right='on', left='on')
    axe.grid()
    pass

runs = range(1, 22)
r2 = range(601, 622)
MSDs = np.empty(shape=(len(runs), 2, 93))
Fs, Cb, V = [np.empty_like(MSDs) for _ in range(3)]
for i, r in enumerate(runs):
    data = np.genfromtxt(f'../T0.13_MC/run_{r}/results/obs.txt', delimiter='').T
    # data2 = np.genfromtxt(f'../IS/T0.12_MC/run_{r}/results/obs.txt', delimiter='').T
    MSDs[i, 0] = data[3]
    V[i, 0] = data[2]
    Fs[i, 0] = data[4]
    Cb[i, 0] = data[5]

    # MSDs[i, 1] = data2[3]
    # V[i, 1] = data2[2]
ts1 = data[0]
print(len(ts1))

MSDs2 = np.empty(shape=(len(r2), 2, 93))
Fs2, Cb2, V2 = [np.empty_like(MSDs2) for _ in range(3)]
for i, r in enumerate(r2):
    data = np.genfromtxt(f'../T0.04_SWAP/run_{r}/results/obs.txt', delimiter='').T
    MSDs2[i, 0] = data[3]
    V2[i, 0] = data[2]
    Fs2[i, 0] = data[4]
    Cb2[i, 0] = data[5]

ts2 = data[0]
print(len(ts2))

plot = sys.argv[1]
if plot == 'u':
    plt.plot(ts1, np.mean(V, axis=0)[0], 'r-', lw=1.2, label='MC')
    plt.plot(ts2, np.mean(V2, axis=0)[0], 'b-', lw=1.2, label='SWAP')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'U/N')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig('T0.13/u.jpg', dpi=300)


elif plot == 'Fs':
    plt.plot(ts1, np.mean(Fs, axis=0)[0], 'r-', lw=1.2, label='MC')
    plt.plot(ts2, np.mean(Fs2, axis=0)[0], 'b-', lw=1.2, label='SWAP')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'$F_s$')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig('T0.13/Fs.jpg', dpi=300)

elif plot == 'Cb':
    plt.plot(ts1, np.mean(Cb, axis=0)[0], 'r-', lw=1.2, label='MC')
    plt.plot(ts2, np.mean(Cb2, axis=0)[0], 'b-', lw=1.2, label='SWAP')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'$C_B$')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig('T0.13/Cb.jpg', dpi=300)

elif plot == 'MSD':
    ax = plt.gca()
    plt.plot(ts1[1:], np.mean(MSDs, axis=0)[0,1:], 'r-', lw=1.2, label='MC')
    plt.plot(ts2[1:], np.mean(MSDs2, axis=0)[0,1:], 'b-', lw=1.2, label='SWAP')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'MSD$(t)$')

    make_axis(ax, True, True)
    plt.legend(frameon=False)
    plt.savefig('T0.13/MSD.jpg', dpi=300)