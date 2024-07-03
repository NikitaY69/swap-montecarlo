import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import sys
import os

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

runs = range(606, 702)
PATH='T0.035'
try:
    os.mkdir(PATH)
except FileExistsError:
    pass

ts = np.unique(np.logspace(0, 8, 200, dtype=int))
MSDs = np.empty(shape=(len(runs), 2, len(ts)))
Fs, Cb, V = [np.empty_like(MSDs) for _ in range(3)]
for i, r in enumerate(runs):
    data = np.genfromtxt(f'../{PATH}/run_{r}/results/obs.txt', delimiter='').T
    # data2 = np.genfromtxt(f'../IS/{PATH}/run_{r}/results/obs.txt', delimiter='').T
    MSDs[i, 0] = data[3]
    V[i, 0] = data[2]
    Fs[i, 0] = data[4]
    Cb[i, 0] = data[5]

    # MSDs[i, 1] = data2[3]
    # V[i, 1] = data2[2]
ts = data[0]

plot = sys.argv[1]
if plot == 'u':
    plt.plot(ts, np.mean(V, axis=0)[0], 'b-', lw=1.2)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'U/N')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig(f'{PATH}/u.jpg', dpi=300)


elif plot == 'Fs':
    plt.plot(ts, np.mean(Fs, axis=0)[0], 'b-', lw=1.2)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'$F_s$')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig(f'{PATH}/Fs.jpg', dpi=300)

elif plot == 'Cb':
    plt.plot(ts, np.mean(Cb, axis=0)[0], 'b-', lw=1.2)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'$C_B$')

    ax = plt.gca()
    make_axis(ax, True, False)
    plt.legend(frameon=False)
    plt.savefig(f'{PATH}/Cb.jpg', dpi=300)

elif plot == 'MSD':
    ax = plt.gca()
    plt.plot(ts[1:], np.mean(MSDs, axis=0)[0,1:], 'r-', lw=1.2, label='Physical trajectory')
    # plt.plot(ts[1:], np.mean(MSDs, axis=0)[1,1:], 'b-', lw=1.2, label='Inherent trajectory')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time step')
    plt.ylabel(r'MSD$(t)$')

    make_axis(ax, True, True)
    plt.legend(frameon=False)
    plt.savefig(f'{PATH}/MSD.jpg', dpi=300)