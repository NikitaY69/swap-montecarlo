from smc_db import RunsFactory
import numpy as np

from scipy.sparse.linalg import eigs
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable

import soft_spot_project.sim_lib as sim
from soft_spot_project.init_model import InitializeModel
from soft_spot_project.modes_tools import get_sparse_hessian,participation_ratio

from softspot.softspot import SoftSpot

class SoftModes(RunsFactory):
    def __init__(self, root, idx, t):
        RunsFactory.__init__(self, root)
        self.run = self.db[idx]
        self.N = self.run['N']
        self.L = np.sqrt(self.N)
        self.T = self.run['T']
        self.set_snapshot(t)
        self.box_size = self.L*np.ones(self.cfg.shape[1])
        sim = InitializeModel(model_name='swap_ipl',r=self.cfg,diameter=self.sigma,box_size=self.box_size,strain=0.0).sim

        # compute "everything" (information needed for the forces and hessian), has to be called before any mode analysis
        sim.pairwise.compute_poly_everything()
        print('###################################################################')
        print('U/N=',sim.system.thermo_pot/sim.system.npart)
        print('Pressure (excess)=',sim.system.thermo_pre)
        print('Stress (xy)=',sim.system.thermo_sigma)
        print('Typical force grad =',sim.system.typical_grad)

    def get_modes(self, n=50, harmonic=False):
        n_modes = n + sim.system.ndim
        # compute Hessian
        self.hessian = get_sparse_hessian()
        # partial diagonalization
        vals, vecs = eigs(self.hessian, k=n, which='LM', sigma=0, tol=0)
        # Throw away infinitesimal imaginary part and sort.
        vals  = np.real(vals)
        vecs = np.real(vecs)
        idx = vals.argsort()
        vals   = vals[idx]
        vecs  = vecs[:,idx]

        softspot = SoftSpot(ndim=sim.system.ndim, npart=sim.system.npart, hessian=self.hessian)

        pi_modes = []
        pi_kappa = []
        ep_psi = []
        ep_pi  = []

        count_pi = 0

        for k in range(2,len(vals)):
            result_cg = softspot.find(vecs[:,k], mode='cg', options={'tol': 1e-6})
            pi = result_cg['pi']
            kappa_pi = result_cg['kappa']
            ep_psi.append(participation_ratio(vecs[:,k]))

            already_found = False
            for m,mode in enumerate(pi_modes):
                dot = np.abs(np.dot(pi,mode))
                if(dot>0.999):
                    already_found=True
                    break

            if (already_found==False):
                pi_modes.append(pi)
                pi_kappa.append(kappa_pi)
                ep_pi.append(participation_ratio(pi))
                count_pi+=1
                print('#',k,'pi count ',count_pi,'kappa(pi)',kappa_pi)


        pi_modes  = np.array(pi_modes)
        pi_kappa  = np.array(pi_kappa)
        psi_kappa = vals[2:]
        ep_psi = np.array(ep_psi)
        ep_pi  = np.array(ep_pi)

        # sort according to stiffness
        idx = pi_kappa.argsort()
        pi_modes = pi_modes[idx]
        pi_kappa = pi_kappa[idx]
        ep_pi = ep_pi[idx]

        if harmonic:
            return psi_kappa, vecs, ep_psi
        else:
            return pi_kappa, pi_modes, ep_pi
        
    def set_snapshot(self, t):
        data = np.loadtxt(f"{self.run['rootdir']}/configs/cfg_{t}.xy")
        self.sigma = data[:,0]
        self.cfg = self.shift_box(data[:,1:])

    def shift_box(self, a):
        return a - self.L * np.floor((a + (self.L/2)) / (self.L)) + self.L/2