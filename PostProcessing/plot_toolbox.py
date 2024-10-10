import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.animation as animation
from functools import partial
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from smc_db import RunsFactory
from colorama import Fore, Style

class PlotToolBox(RunsFactory):
    def __init__(self, root, idx=False, figsize=(7,7), fontsize=20, linewidth=2.5, facecolor='black', cbar=True):
        '''
        This class gives tools to plot various quantities (observables aside) 
        of interest in the context of SMC runs.
        '''
        RunsFactory.__init__(self, root)
        if idx == False:
            print(f'{Fore.YELLOW}You have not provided an index for the database to point at...{Style.RESET_ALL}\n' + 
                  f'Database by default pointing at idx {Fore.RED}0{Style.RESET_ALL}.\n' + 
                  f'Use {Fore.CYAN}self.show(){Style.RESET_ALL} to get all database details.\n')
            idx = 0
        self.run = self.db[idx]
        self.N = self.run["N"]
        self.L = np.sqrt(self.N)/2
        self.lin_ts = np.linspace(0, self.run['steps'], self.run['linPoints'], endpoint=False, dtype=int)
        self.lin_ts[0] += 1
        self.log_ts = np.unique(np.logspace(0, np.log10(self.run['steps']), self.run['logPoints'], dtype=int))
        self.fig, self.ax, self.cbar_ax = self.make_canvas(figsize, fontsize, linewidth, facecolor, cbar)
        self.multi_render = False

    @staticmethod
    def make_canvas(figsize=(7,7), fontsize=20, linewidth=2.5, facecolor='black', cbar=True):
        # Plot params
        plt.rcParams['figure.figsize'] = figsize
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.linewidth'] = linewidth
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['animation.embed_limit'] = 2**128

        # Plot setup
        fig, ax = plt.subplots(dpi=500)
        ax.set_facecolor(facecolor)
        cbar_ax = None
        if cbar:
            divider = make_axes_locatable(ax)
            cbar_ax = divider.append_axes("right", size="5%", pad=0.1)

        return fig, ax, cbar_ax

    def init_fig(self, particles=False, disp_field=False, \
                 c_p='grey', alpha_p=0.4, lw_p=0.5, \
                 cmap_f='jet', scale=1, width=0.00375/1.5, headwidth=2.5, headlength=2.5, headaxislength=2.5):
        # Ici plutôt créer des dictionnaires de keywords pour chaque champ
        if not particles and not disp_field:
            raise NotImplementedError(f'{Fore.YELLOW}You are asking me to do nothing...{Style.RESET_ALL}\n' + 
                                      f'You need to set at least one plot, ' + 
                                      f'either {Fore.CYAN}pos{Style.RESET_ALL} or {Fore.CYAN}field{Style.RESET_ALL}.')
        if particles:
            self.c0 = self.load_cfg(t=1)
            x, y, r = self.c0; x_box, y_box = [self.Pshift(a) for a in [x, y]]
            dX, dY = [[0]*self.N for _ in range(2)]
            if self.cbar_ax == None:
                self.particles = [Circle((x_box[i], y_box[i]), r[i], facecolor=c_p, linewidth=lw_p, alpha=alpha_p) for i in range(self.N)]
                self.collection = PatchCollection(self.particles, match_original=True)
            else:
                self.particles = [Circle((x_box[i], y_box[i]), r[i], linewidth=lw_p, alpha=alpha_p) for i in range(self.N)]
                self.collection = PatchCollection(self.particles, cmap=c_p, match_original=True)
            self.ax.add_collection(self.collection)
            # scat = self.ax.scatter([], [], s=[], color=c_p, alpha=alpha_p, linewidths=lw_p)
        if disp_field:
            self.quiv = self.ax.quiver(x_box, y_box, # XY\ 
                                  dX, dY, # UV\
                                  np.hypot(dX, dY), # C\
                                  cmap=cmap_f, scale=scale, scale_units='xy', angles='xy', \
                                  width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
        self.ax.set_xlim([-self.L, self.L])
        self.ax.set_ylim([-self.L, self.L])

    def render_stuff(self, t, particles=False, disp_field=False, A_p=None, **kwargs):
        if not self.multi_render:
            self.init_fig(particles, disp_field, **kwargs)
        if t != 1:
            self.ax.set_title(f"t={t}")
        else:
            self.ax.set_title(f"t=0")
        x, y, r = self.load_cfg(t)
        x_box, y_box = [self.Pshift(a) for a in [x,y]]
        if particles:
            for i, p in enumerate(self.particles):
                p.set_center((x_box[i], y_box[i]))
                p.set_radius(r[i])
            self.collection.set_paths(self.particles)
            if A_p != None:
                ind_t = np.argwhere(self.lin_ts == t)
                self.collection.set_array(A_p[ind_t])
                self.collection.set_clim(A_p[ind_t].min(), A_p[ind_t].max())
        if disp_field:
            dX = x - self.c0[0]
            dY = y - self.c0[1]
            dR = np.hypot(dX, dY)
            # xy = np.vstack((x_box, y_box)).T
            # self.quiv.set_offsets(xy)
            self.quiv.set_UVC(dX, dY, np.hypot(dX, dY))
            self.quiv.norm.vmin = dR.min()
            self.quiv.norm.vmax = dR.max()
        if not self.multi_render:
            plt.savefig("test2.jpg")
        else:
            return self.collection, self.quiv,

    def movie(self, save, log=False, n_frames = 'all', fps=8, interval=50, **kwargs):
        self.multi_render = True
        self.init_fig(**kwargs)
        if not log:
            ts = self.lin_ts
        else:
            ts = self.log_ts
        if n_frames != 'all':
            ts = ts[:n_frames]
        flow = animation.FuncAnimation(self.fig, partial(self.render_stuff, **kwargs), frames=ts, interval=interval, blit=True)
        writermp4 = animation.FFMpegWriter(fps=fps)
        flow.save(f"{save}", writer=writermp4)

    def Pshift(self, a):
        # Periodic shift in main box of size L
        return a - 2 * self.L * np.floor((a + self.L) / (2*self.L))

    def load_cfg(self, t):
        t = int(t)
        diameter, x, y = np.genfromtxt(f"{self.run['rootdir']}/configs/cfg_{t}.xy").T
        x -= np.mean(x); y -= np.mean(y) # Substracting center of mass movement
        return x, y, diameter/2