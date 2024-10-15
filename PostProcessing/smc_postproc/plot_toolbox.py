# Numpy
import numpy as np
# Matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.animation as animation
from functools import partial
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
# Runs
from smc_db import RunsFactory
# Colorama
from colorama import Fore, Style

class PlotToolBox(RunsFactory):
    """
    This class gives tools to plot various quantities (observables aside) 
    of interest in the context of SMC runs. For example, by just loading 
    a specific run, you can automatically get the movie of its evolution, 
    render a specific snapshot, etc.
    """
    def __init__(self, root, idx=False, figsize=(7,7), facecolor='black', cbar=True, \
                 cbar_lim=None, cbar_label=None):
        '''
        Initialization of canvas object.

        Parameters
        ----------
        root : str
            Path to RunsFactory pkl root

        idx : int
            Index of run in consideration

        figsize : tuple
            Figsize
        
        facecolor : float
            Color of the canvas background

        cbar : bool
            whether or not to add a colorbar
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
        self.lin_ts = np.linspace(0, self.run['tau'], self.run['linPoints'], endpoint=False, dtype=int)
        self.lin_ts[0] += 1
        self.log_ts = np.unique(np.logspace(0, np.log10(self.run['tau']), self.run['logPoints'], dtype=int))

        # Fig attributes
        self.fig, self.ax, self.cbar_ax = self.make_canvas(figsize=figsize, facecolor=facecolor, cbar=cbar)
        self.cbar_lim = cbar_lim
        self.cbar_label = cbar_label
        self.multi_render = False

    @staticmethod
    def make_canvas(figsize=(7,7), fontsize=15, linewidth=2.5, facecolor='black', cbar=True):
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

    def init_fig(self, particles, disp_field, A_p, plot_params):
        """
        Init. fig specifities.

        Parameters
        ----------
        particles : bool
            whether or not to plot particles

        disp_field : bool
            whether or not to plot displacement fields

        A_p : array of shape (self.ts, self.N)
            Mapping array for particles color (by default None)

        plot_params : dic
            Dictionnary of plotting params (see just below for infos)
        """
        # Particle-related params
        c_p = plot_params.get('c_p', 'grey')            # uniform color
        cmap_p = plot_params.get('cmap_p', 'Blues')     # if not, cmap
        alpha_p = plot_params.get('alpha_p', 0.4)       # opacity
        lw_p = plot_params.get('lw_p', 0.5)             # circles' linewidth
        
        # Field-related params
        cmap_f = plot_params.get('cmap_f', 'jet')       # field norm cmap 
        scale = plot_params.get('scale', 1)             # scaling with axis units
        # Arrow specs
        width = plot_params.get('width', 0.00375 / 1.5)
        headwidth = plot_params.get('headwidth', 2.5)
        headlength = plot_params.get('headlength', 2.5)
        headaxislength = plot_params.get('headaxislength', 2.5)

        if not particles and not disp_field:
            raise NotImplementedError(f'{Fore.YELLOW}You are asking me to do nothing...{Style.RESET_ALL}\n' + 
                                      f'You need to set at least one plot, ' + 
                                      f'either {Fore.CYAN}pos{Style.RESET_ALL} or {Fore.CYAN}field{Style.RESET_ALL}.')
        if particles:
            self.c0 = self.load_cfg(t=1) # init configuration
            x, y, r = self.c0; x_box, y_box = [self.Pshift(a) for a in [x, y]]
            dX, dY = [[0]*self.N for _ in range(2)]
            
            if A_p is None:
                self.particles = [Circle((x_box[i], y_box[i]), r[i], \
                                         facecolor=c_p, linewidth=lw_p, alpha=alpha_p)
                                  for i in range(self.N)]
                self.collection = PatchCollection(self.particles, match_original=True)
            else:
                self.particles = [Circle((x_box[i], y_box[i]), r[i], \
                                         linewidth=lw_p, alpha=alpha_p) 
                                  for i in range(self.N)]
                self.collection = PatchCollection(self.particles, cmap=cmap_p, match_original=True)
                cbar = self.fig.colorbar(self.collection, cax=self.cbar_ax, label=self.cbar_label)

            self.ax.add_collection(self.collection)

        if disp_field:
            self.quiv = self.ax.quiver(x_box, y_box, # XY
                                dX, dY,              # UV
                                np.hypot(dX, dY),    # C
                                cmap=cmap_f, scale=scale, scale_units='xy', angles='xy', \
                                width=width, headwidth=headwidth, \
                                headlength=headlength, headaxislength=headaxislength)
            
        self.ax.set_xlim([-self.L, self.L])
        self.ax.set_ylim([-self.L, self.L])

    def render_stuff(self, t, particles=False, disp_field=False, A_p=None, plot_params={}, save=None):
        """
        Render snapshots with displacement fields.
        This function can be used both for single-rendering or multi-rendering
        (as an update function for FuncAnimation).

        Parameters
        ----------
        t : int
            Saved timestep for run at hand.

        particles, disp_field, A_p, plot_params : 
            check self.init_fig

        save : str
            Path to saving file
        """
        if not self.multi_render:
            self.init_fig(particles, disp_field, A_p, plot_params)

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

            if A_p is not None:
                ind_t = np.argwhere(self.lin_ts == t)[0]
                self.collection.set_array(np.squeeze(A_p[ind_t]))
                if self.cbar_lim is None:
                    c_min, c_max = A_p[ind_t].min(), A_p[ind_t].max()
                else:
                    c_min, c_max = self.cbar_lim
                self.collection.set_clim(c_min, c_max)
                
        if disp_field:
            dX = x - self.c0[0]
            dY = y - self.c0[1]
            dR = np.hypot(dX, dY)
            self.quiv.set_UVC(dX, dY, np.hypot(dX, dY))
            self.quiv.norm.vmin = dR.min()
            self.quiv.norm.vmax = dR.max()

        if not self.multi_render:
            if save is None:
                plt.show()
            else:
                plt.savefig(f"{save}")

        else:
            # In case used as an update function for FuncAnimation
            return self.collection, self.quiv,

    def movie(self, particles=False, disp_field=False, A_p=None, \
              log=False, n_frames = 'all', fps=8, \
              plot_params={}, save=None):
        """
        Makes an animation of a run.

        Parameters
        ----------
        particles, disp_field, A_p: 
            check init_fig

        log : bool
            whether or not to use lin-spaced or log-spaced saved configurations
        
        n_frames = int
            number of snapshots to be considered (by default 'all')

        fps = int
            Frames per second

        plot_params : dic
            check init_fig

        save : str
            Path to saving file
        """
        self.multi_render = True
        self.init_fig(particles, disp_field, A_p, plot_params)

        if not log:
            ts = self.lin_ts
        else:
            ts = self.log_ts 
        if n_frames != 'all':
            ts = ts[:n_frames]

        flow = animation.FuncAnimation(self.fig, partial(self.render_stuff, 
                                                         particles=particles, 
                                                         disp_field=disp_field, 
                                                         A_p=A_p, 
                                                         plot_params=plot_params), 
                                        frames=ts, interval=50, blit=True)
        writermp4 = animation.FFMpegWriter(fps=fps)

        if save == None:
            plt.show()
        else:
            flow.save(f"{save}", writer=writermp4)

    def Pshift(self, a):
        # Periodic shift in main box of size L
        return a - 2 * self.L * np.floor((a + self.L) / (2*self.L))

    def load_cfg(self, t):
        """
        Load cfg file for a specific timestep.
        """
        t = int(t)
        diameter, x, y = np.genfromtxt(f"{self.run['rootdir']}/configs/cfg_{t}.xy").T
        x -= np.mean(x); y -= np.mean(y) # Substracting center of mass movement
        return x, y, diameter/2