from smc_db import RunsFactory
import numpy as np

class Observables(RunsFactory):
    def __init__(self, root):
        '''
        This class computes ensemble observables given a database
        '''
        RunsFactory.__init__(self, root)
        # self.get_ensemble(self, params)
        self.set_run(self.ensemble[0]) # setting first run to define some params
        self.lin_ts = np.linspace(0, self.run['steps'], self.run['linPoints'], endpoint=False)
        self.lin_ts[0] += 1
        self.log_ts = np.unique(np.logspace(0, self.run['steps'], self.run['logPoints']), dtype=int)

    def compute_average(self, obs, log=True):
        if log:
            self.ts = self.lin_ts
        else:
            self.ts = self.log_ts
        obs_ = np.empty(shape=(len(self.ensemble), len(self.ts)))
        for i, run in enumerate(self.ensemble):
            obs_[i] = np.genfromtxt(f'{run['rootdir']}/obs.txt', delimiter='', \
                                    usecols=[run[obs]])
        
        return np.mean(obs_, axis=0)

    def set_ensemble(self, idx):
        self.ensemble = self.db[idx]

    def get_ensemble(self, params):
        return None