from smc_db import RunsFactory
from modes import SoftModes
# from observables import Observables
from soft_spot_project.visualization_tools import render_snapshot,render_field
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.use('Agg')

import numpy as np

db = RunsFactory('test.pkl')
# db.show()
db.insert('/home/allaglo/production/IS/T0.04_SWAP/run_601/results')
db.set_run(0)
db.show()
# obs = Observables('test.pkl')
# obs.set_ensemble([0])
# print(obs.lin_ts)
# MSD = obs.compute_average('MSD')
# print(MSD)

# PH_modes = SoftModes('test.pkl', 0, 1)
# freq, modes, PR = PH_modes.get_modes(n=10)
# ax = plt.gca()
# render_snapshot(ax,PH_modes.cfg,PH_modes.sigma,PH_modes.box_size,color='dimgray',edgecolor='black',alpha=0.75)

# for k,mode in enumerate(modes): # display the softest at the end
# 	render_field(ax,PH_modes.cfg,mode,scale=0.02, remove_fraction=0.95)
	
# plt.savefig('test.jpg', dpi=500)