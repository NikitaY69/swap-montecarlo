from .smc_db import RunsFactory
import numpy as np
from colorama import Fore, Style

class Observables(RunsFactory):
    def __init__(self, root, ensemble='all'):
        '''
        This class computes ensemble observables given a database
        '''
        RunsFactory.__init__(self, root)
        self.set_ensemble(ensemble)

    def compute_average(self, obs, log=True):
        if log:
            self.ts = self.log_ts
        else:
            self.ts = self.lin_ts
        obs_ = np.empty(shape=(len(self.ensemble), len(self.ts)))
        for i, run in enumerate(self.ensemble):
            obs_[i] = np.genfromtxt(f'{run["rootdir"]}/obs.txt', delimiter='', \
                                    names=True)[obs]
        
        return np.mean(obs_, axis=0)

    def set_ensemble(self, idx='all'):
        '''
        Defining ensemble over which observables are averaged.
        '''
        self.ensemble = np.array(self.db, dtype=dict)
        if type(idx) is int:
            raise IndexError(f'When retrieving single-run observables, ' +
                             f'you must still provide the index in a list, '+
                             f'ie {Fore.CYAN}[0]{Style.RESET_ALL} instead of {Fore.RED}0{Style.RESET_ALL}.')
        if idx != 'all':
            self.ensemble = self.ensemble[idx]
        dummy = self.db[0]
        
        # Timesteps 
        self.lin_ts = np.linspace(0, dummy['tau'], dummy['linPoints'], endpoint=False)
        self.lin_ts[0] += 1
        self.log_ts = np.unique(np.logspace(0, np.log10(dummy['tau']), dummy['logPoints'], dtype=int))

    # def get_ensemble(self, params):
    #     return None