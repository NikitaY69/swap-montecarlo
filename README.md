# Swap Monte Carlo algorithm for 2-dimensional polydisperse liquids

This repository provides a C++ implementation of Swap Monte Carlo (SMC) simulations
for 2-dimensional continuous polydisperse liquids. For informations about the model and the algorithm used, please refer to [A. Ninarello et al., Phys. Rev. X 7, 021039 (2019)](https://link.aps.org/doi/10.1103/PhysRevX.7.021039) and [L. Berthier et al., J. Stat. Mech. (2019) 064004](https://iopscience.iop.org/article/10.1088/1742-5468/ab1910). 

A python post-processing package is also made available. One can refer to `PostProcessing/README.md` for instructions for using the later.

## How to use it ?
Our SMC module is usable for different purposes as can testify the various branches avalailable. The `master` branch can be used for both thermalization and production runs. The `observables-only` branch is made to compute observables AFTER a simulation has finished. A `continue` branch finally manages to carry on finished runs. On each branch, the module can be compiled with
```
cd src/SMC/
g++ *.cpp -o EXEC_NAME -lstdc++fs -O2 -mfma -mbmi2 -flto -lboost_program_options
```

P.S. once you have compiled an executable under each relevent branch, you no longer need to switch between branches as you can just continuously re-use your executables.

### Thermalization
Thermalization is ensuring the fact that starting from any initial configuration, equilibrium is reached at desired temperature. One can run thermalization runs with
```
EXEC_NAME --input INPUT_FILE \ 
          --outdir OUT_DIRECTORY \
          --N 5 \
          --T 0.04
          --tau 100000 \
          --tw 1 \
          --cycles 1 \
          --lin 50 \
          --log 50 \
          --p_swap 0.2 \
          [--MSD --Cb --Fs --U]
```

- `input`: path to starting configuration with data structure `X Y SIGMA` repeated `N` times (see `tutorial/cfg.xy` for example)
- `outdir`: path of output directory
- `N`: number of interacting particles
- `T`: temperature in reduced units
- `tau`: single-run time (as in aging plots, time after which observables are reset)
- `tw`: waiting-time between each cycle (as in aging plots, time after which a new cycle is measured)
- `cycles`: number of cycles of decorrelation tours
- `lin`: number of lin-spaced checkpoints (configurations)
- `log`: number of log-spaced checkpoints (configurations AND observables)
- `p_swap`: SWAP probability
- `MSD`, `Cb`, `Fs`, `U`: flags to compute the _Mean-Squared Displacement_, the _bond-breaking correlation function_, the _self-part of the intermediate scattering function_ and the _average potential energy_


The total number of steps is calculated as in $n_\texttt{steps}=\texttt{tw}*(\texttt{cycles}-1)+\tau$.
To check-on equilibration, one must use $\texttt{tw}>1$ and $\texttt{cycles}>1$. Usually, one launches a first run with a `tw`=`cycles`=1 to have an approximate of the relaxation time $\tau_\alpha$ with $C_B$ or $F_s$. A good value for `tw` is then slightly below $\tau_\alpha$ using `cycles`>1; finally one reaches equilibrium when 1- no memory-effects are visible between different cycles and 2- the potential energy is constant.

### Production
Calculating most physical observables relies on having good ensemble statistics. To do so, after reaching equilibrium at desired temperature, one must get observables evolution from different starting configurations. To proceed, you can use the exact same command as in [Thermalization](#thermalization) but with $\texttt{tw}=\texttt{cycles}=1$.<br>
Most of the time, you can also use a slurm job with arrays each one running a specific starting configuration.

### Inherent structures
### Observables-only
### Runs continuation

## Workflow example

In `tutorial/tutorial.ipynb`, we provide an example of a complete typical workflow:

1. Run one SMC thermalization run
1. Ensure equilibration
1. Produce multiple equilibrated configurations
1. Visualize the dynamics
1. Compute inherent structures
1. Compute and visualize soft modes