import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
runs = range(601, 702)
PATH = 'T0.035'
# Constants
N = 2000
L = np.sqrt(N) / 2
ss = np.arange(1, 100000000, 500000)
# rs = [35 * k / 49 for k in range(2, 51)]
# i_r = 2
# print(rs[i_r])
i_run = int(sys.argv[1])

# Helper functions
def Pshift(a):
    return a - 2 * L * np.floor((a + L) / (2*L))

# Data loading
cfgs = np.empty(shape=(len(ss), 3, N))
# micros = np.empty(shape=(len(ss), 51, N))
dR = np.empty(shape=(len(ss), 2, N))

for i, t in enumerate(ss):
    if i != 0:
        t-=1
    cfgs[i] = np.genfromtxt(f'../{PATH}/run_{i_run}/results/configs/cfg_{t}.xy').T
    if i > 0:
        cfgs[i, 1] -= np.mean(cfgs[i, 1])
        cfgs[i, 2] -= np.mean(cfgs[i, 2])
    # if i != 0:
    #     micros[i] = np.genfromtxt(f'T0.04/run_{i_run}/results/micro_p/products_loc_{t}.txt').T
    dR[i, 0] = cfgs[i,1]-cfgs[0,1]
    dR[i, 1] = cfgs[i,2]-cfgs[0,2]

dR[0] = np.zeros_like(dR[1])
# micros[0] = np.zeros_like(micros[1])

# Plot setup
fig, ax = plt.subplots(dpi=500)
ax.set_facecolor('black')

# divider = make_axes_locatable(ax)
# cbar_ax = divider.append_axes("right", size="5%", pad=0.1)
vmax = cfgs[0, 0].max()-cfgs[0, 0].min()
scat = ax.scatter([], [], s=[], color='grey', alpha=0.4)
quiver = ax.quiver(np.empty(shape=(N)), np.empty(shape=(N)), np.empty(shape=(N)), np.empty(shape=(N)), \
                   np.empty(shape=(N)), cmap='jet', scale=1, scale_units='xy', angles='xy', width=0.0035, \
                   headwidth=3, headlength=3, headaxislength=2.5)  # 0.00375
pos0 = np.vstack((Pshift(cfgs[0, 1]), Pshift(cfgs[0, 2]))).T
quiver.set_offsets(pos0)

def init():
    ax.set_xlim([-L, L])
    ax.set_ylim([-L, L])
    scat.set_offsets(np.empty((N, 2)))
    scat.set_sizes(np.empty((N)))
    # scat.set_array(np.empty((N)))
    quiver.set_UVC(np.zeros((N)), np.zeros((N)), C=np.zeros(N))

    return scat, quiver,

def update(i):
    print(i, flush=True)
    ax.set_title(f"t={ss[i]-1}")
    # Normalization
    # if i > 0:
    #     divnorm = TwoSlopeNorm(vcenter=0, vmax=micros[i, i_r, :].max(), vmin=micros[i, i_r, :].min())
    #     scat.set_norm(divnorm)
    pos = np.vstack((Pshift(cfgs[i, 1]), Pshift(cfgs[i, 2]))).T
    scat.set_offsets(pos)
    scat.set_sizes(60 * cfgs[i, 0] ** 2)
    # scat.set_array(micros[i, i_r, :])
    # cbar = fig.colorbar(scat, cax=cbar_ax, label=r'$\langle\Delta\mathbf{r}_i\cdot\Delta\mathbf{r}_j\rangle_{j-i\leq r}$')

    C = np.hypot(dR[i, 0], dR[i, 1])
    quiver.set_UVC(dR[i, 0], dR[i, 1], C)
    quiver.norm.vmin = C.min()
    quiver.norm.vmax = C.max()

    return scat, quiver

anim = animation.FuncAnimation(fig, update, init_func=init, frames=len(ss), interval=50, blit=True)
writermp4 = animation.FFMpegWriter(fps=3)
anim.save(f"../{PATH}/run_{i_run}/disp.mp4", writer=writermp4)