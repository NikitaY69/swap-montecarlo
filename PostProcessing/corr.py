import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import matplotlib.ticker as mticker
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
PATH = 'T0.035'
try:
    os.mkdir(PATH)
except FileExistsError:
    pass

N = 2000
L = np.sqrt(N)/2
sigma_m = 1.000218223
sigmaMax = 1.613048
ts = np.unique(np.logspace(0, 8, 200, dtype=int))
zs = np.empty(shape=(len(runs), len(ts)))
rs = [35*k/49 for k in range(2,50)]

def Pshift(a):
    return a-2*L*np.floor((a+L)/(2*L))

for j, r in enumerate(runs):
    if j == 13:
        pass
    else:
        print(j)
        corr = np.genfromtxt(f'../{PATH}/run_{r}/results/products.txt')
        ts = corr.T[0].astype(int)
        MSD = np.genfromtxt(f'../{PATH}/run_{r}/results/obs.txt').T[3]
        # cfgs = {t:None for t in corr.T[0].astype(int)}
        # cfg0 = np.genfromtxt(f'../{PATH}/run_{r}/results/configs/cfg_1.xy').T

        # Load data
        # for t in corr.T[0].astype(int):
        #     cfgs[t] = np.genfromtxt(f'../{PATH}/run_{r}/results/configs/cfg_{t:d}.xy').T
        for i, ss in enumerate(ts):
            if i>6:

                # MSD
                # dX = cfgs[ss][1]-cfg0[1]
                # dY = cfgs[ss][2]-cfg0[2]
                # dX -= np.mean(dX)
                # dY -= np.mean(dY)
                # MSD = np.mean(dX*dX+dY*dY)
                # print(MSD)

                plt.clf()
                c = corr[np.argwhere(corr==ss)[0,0]][3:]/N
                # print(c/MSD[i])
                f = (c-c[-1])/(c[0]-c[-1])
                interp = interp1d(np.array(rs)/(sigma_m), c/MSD[i], kind='cubic')
                def f_(x):
                    return interp(x)-(c.max()/MSD[i])*np.exp(-1)
                zeta = brentq(f_, 2, 25)
                zs[j, i] = zeta

                # ax = plt.gca()
                # xmin, xmax = ax.get_xlim()
                # ymin, ymax = ax.get_ylim()
                # plt.xlim(0, xmax)
                # plt.ylim(ymin,ymax)
                # plt.vlines(zeta, ymin=ymin, ymax=(c.max()/MSD[i])*np.exp(-1), colors='k', ls='--', lw=2)
                # plt.hlines((c.max()/MSD[i])*np.exp(-1), xmin=0, xmax=zeta, colors='k', ls='--', lw=2)
                # make_axis(ax, False)
                # plt.grid(True)
                # plt.xlabel(r'$\frac{r}{\bar{\sigma}}$')
                # norm = r'$\frac{f(r,t)-f(r_\mathrm{max},t)}{f(r_\mathrm{min},t)-f(r_\mathrm{max},t)}$'
                # plt.ylabel(r'$\frac{f(r,t)}{\mathrm{MSD}}$')

                # plt.text(20, 0.3, r'$f(r,t)=\langle\Delta\mathbf{r}_i\cdot\Delta\mathbf{r}_j\rangle_{j-i\leq r}$')
                # plt.text(4, 0.3, r'$\mathrm{e}^{-1}$')
                # plt.savefig('figs/corrs/corr_%d.jpg' %ss, dpi=300)
        
plt.cla()
plt.clf()
ax = plt.gca()

plt.xlabel(r'Time step ($\times 10^6$)')
plt.ylabel(r'$\zeta$')
plt.plot(ts[7:], np.mean(zs, axis=0)[7:len(ts)], 'b-', lw=1.2)
plt.xscale('log')
plt.yscale('log')
formatter = mticker.ScalarFormatter()
ax.yaxis.set_minor_formatter(formatter)
# ax.ticklabel_format(style='plain', axis='y')
make_axis(ax, True, True)
plt.grid(True)
plt.grid(True, which='minor', axis='y')
plt.savefig(f'{PATH}/zeta.jpg', dpi=300)
