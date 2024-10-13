# Swap Monte Carlo algorithm for 2-dimensional polydisperse liquids

This repository provides a C++ implementation of Swap Monte Carlo (SMC) simulations
for 2-dimensional continuous polydisperse liquids with unit density. The SMC algorithm consists in adding a supplementary trial move: exchanging particle diameters. The latter has proven to drastically speed up glassy dynamics when combined with casual displacement moves. <br>
For informations about the model and the algorithm used, please refer to [A. Ninarello _et al._, Phys. Rev. X 7, 021039 (2019)](https://link.aps.org/doi/10.1103/PhysRevX.7.021039) and [L. Berthier _et al._, J. Stat. Mech. (2019) 064004](https://iopscience.iop.org/article/10.1088/1742-5468/ab1910). 

A python post-processing package is also made available: one can refer to `PostProcessing/README.md` for instructions for using it.

## Compilation and execution
On each branch, the module can be compiled with
```
cd src/SMC/
g++ *.cpp -o EXEC_NAME -lstdc++fs -O2 -mfma -mbmi2 -flto -lboost_program_options
```
The executable is then parsed as follows
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


Because the model has a density of one, the size of the system is simply taken as $L=\sqrt{N}$. Particles then interact with each other in a box $[-L/2,L/2]\times [-L/2,L/2]$ with the minimum image convention. <br>
The chosen unit of time is one monte-carlo **sweep** (that is $N$ consecutive displacement/swap trials). The total number of steps is calculated as in $n_\texttt{steps}=\texttt{tw}*(\texttt{cycles}-1)+\tau$.

## How to use it ?
Our SMC module is usable for different purposes as can testify the various branches avalailable. The `master` branch can be used for both thermalization and production runs. The `observables-only` branch is made to compute observables AFTER a simulation has finished. A `continue` branch finally manages to carry on finished runs. 

P.S. once you have compiled an executable under each relevent branch, you no longer need to switch between branches as you can just continuously re-use your executables.

### Thermalization
Thermalization is ensuring the fact that starting from any initial configuration, equilibrium is reached at desired temperature. <br>
To check on equilibration, one must use $\texttt{tw}>1$ and $\texttt{cycles}>1$. Usually, one launches a first run with a `tw`=`cycles`=1 to have an approximate of the relaxation time $\tau_\alpha$ with $C_B$ or $F_s$. A good value for `tw` is then slightly below $\tau_\alpha$ using `cycles`>1; finally one reaches equilibrium when 1- no memory-effects are visible between different cycles and 2- the potential energy is constant.

### Production
Calculating most physical observables relies on having good ensemble statistics. To do so, after reaching equilibrium at desired temperature, one must get observables evolution from different starting configurations. To proceed, you must use $\texttt{tw}=\texttt{cycles}=1$ and select the observables of interest.<br>
Most of the time, you can also use a slurm job with arrays each one running a specific starting configuration.

### Inherent structures
Conjuguate gradient algorithms are standard to obtain inherent structures (local minima in the potential energy). In `src/IS/`, we provided the source-code of such algorithm adapted to our 2-dimensional model. <br>
The code is compiled with
```
cd src/IS/
gcc *.c -o EXEC_NAME -lm
```
and then executable with
```
EXEC_NAME INPUT_CFG OUTPUT_CFG
```
- `INPUT_CFG`: input configuration (same format as input in [SMC](#compilation-and-execution))
- `OUTPUT_CFG`: output file to write inherent structure

By contrast to [SMC](#compilation-and-execution), you will notice that no argument parser is provided here - the reason being that the code was not written by me, but by Andrea Ninarello and I only modified it to respect our model. 

A downside for this is that **you must manually enter the number of particles and size of the system in the source-code**. Go to line 25-26 and you will see
```
#define SIZE 20000
#define SIDE 100.0
```
- $\texttt{SIZE}=N\cdot d$ [$d$ being the number of dimensions]
- $\texttt{SIDE}=\sqrt{N}$ with our model

### Observables-only

### Runs continuation
_This branch is still under development. Do not use it yet: an update will come soon._

## Workflow example

In `tutorial/tutorial.ipynb`, we provide an example of a complete typical workflow:

1. Run one SMC thermalization run
1. Ensure equilibration
1. Produce multiple equilibrated configurations
1. Visualize the dynamics
1. Compute inherent structures
1. Compute and visualize soft modes