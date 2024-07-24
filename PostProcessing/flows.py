import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.animation as animation
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

# Plot params
mpl.rcParams['figure.figsize'] = (7, 7)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.linewidth'] = 2.5
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['animation.embed_limit'] = 2**128

# ARRAY PARAMS
PATH = 'count/T0.025_N10000/run_81'
# Constants
N = 10000
L = np.sqrt(N) / 2
ss = np.arange(1, 10000000, 50000)

# Helper functions
def Pshift(a):
    return a - 2 * L * np.floor((a + L) / (2*L))

# Data loading
cfgs = np.empty(shape=(len(ss), 4, N))
dR = np.empty(shape=(len(ss), 2, N))

for i, t in enumerate(ss):
    if i != 0:
        t-=1
    cfgs[i] = np.genfromtxt(f'../{PATH}/results/configs/cfg_{t}.xy').T
    if i > 0:
        cfgs[i, 1] -= np.mean(cfgs[i, 1])
        cfgs[i, 2] -= np.mean(cfgs[i, 2])
    dR[i, 0] = cfgs[i,1]-cfgs[0,1]
    dR[i, 1] = cfgs[i,2]-cfgs[0,2]

dR[0] = np.zeros_like(dR[1])

# Plot setup
fig, ax = plt.subplots(dpi=500)
ax.set_facecolor('black')
divider = make_axes_locatable(ax)
cbar_ax = divider.append_axes("right", size="5%", pad=0.1)

# i = 50
# scat = plt.scatter(Pshift(cfgs[i, 1]), Pshift(cfgs[i, 2]), s=12*cfgs[i, 0]**2, color='grey', linewidths=0.5, alpha=0.4)
# C = np.hypot(dR[i, 0], dR[i, 1])

# disp = plt.quiver(Pshift(cfgs[0, 1]),Pshift(cfgs[0, 2]), 2*dR[i, 0], 2*dR[i, 1], C, cmap='jet', scale = 1, scale_units = 'xy', angles='xy',\
#                     width=0.00375/1.5, headwidth=2.5, headlength=2.5, headaxislength=2.5)
# plt.savefig('test.jpg', dpi=500)

# divider = make_axes_locatable(ax)
# cbar_ax = divider.append_axes("right", size="5%", pad=0.1)

scat = ax.scatter([], [], s=[], c=[], cmap='gist_gray', alpha=0.7, linewidths=0.5)
quiver = ax.quiver(np.empty(shape=(N)), np.empty(shape=(N)), np.empty(shape=(N)), np.empty(shape=(N)), \
                   np.empty(shape=(N)), cmap='jet', scale=1, scale_units='xy', angles='xy', width=0.00375/1.5, \
                   headwidth=2.5, headlength=2.5, headaxislength=2.5)  # 3 3 2.5 # 0.00375
pos0 = np.vstack((Pshift(cfgs[0, 1]), Pshift(cfgs[0, 2]))).T
quiver.set_offsets(pos0)
cbar = fig.colorbar(scat, cax=cbar_ax, label=r'$n_\mathrm{SWAP}$')
def init():
    ax.set_xlim([-L, L])
    ax.set_ylim([-L, L])
    scat.set_offsets(np.empty((N, 2)))
    scat.set_sizes(np.empty((N)))
    scat.set_array(np.empty((N)))
    quiver.set_UVC(np.zeros((N)), np.zeros((N)), C=np.zeros(N))

    return scat, quiver,

def update(i):
    print(i, flush=True)
    ax.set_title(f"t={ss[i]-1}")    
    pos = np.vstack((Pshift(cfgs[i, 1]), Pshift(cfgs[i, 2]))).T
    scat.set_offsets(pos)
    scat.set_sizes(12 * cfgs[i, 0] ** 2) # 60-12
    if i>0:
        A = cfgs[i, -1]-cfgs[i-1,-1]#/(N*(ss[i]-1))
        scat.set_array(A)
        scat.set_clim(vmin=A.min(), vmax=A.max())

    C = np.hypot(dR[i, 0], dR[i, 1])
    quiver.set_UVC(dR[i, 0], dR[i, 1], C)
    quiver.norm.vmin = C.min()
    quiver.norm.vmax = C.max()

    return scat, quiver

anim = animation.FuncAnimation(fig, update, init_func=init, frames=len(ss), interval=50, blit=True)
writermp4 = animation.FFMpegWriter(fps=8)
anim.save(f"../{PATH}/disp.mp4", writer=writermp4)