# Post-Processing python package for SMC software
This sub-directory of the SMC software provides various tools to analyze finished runs. <br>
It is grounded on a rustic database so that the user can analyze his  simulations with a global perspective.

## Installation
Because the package requires a specific version of Python, we suggest to dedicate a specific virtual environment for its usage
```
conda create -n postproc python=3.6.13
conda activate postproc
pip install --upgrade pip
pip install -r requirements.txt
```
We recommand to install the module with the edit flag so that it automatically tracks remote updates:
```
git clone git@github.com:NikitaY69/swap_montecarlo.git
cd PostProcessing/
pip install -e .
```

## Usage
`smc_postproc` is composed of three classes which deal with database management, loading observables and plotting purposes. <br>
This section only suggests a concise presentation of their functionalities: the interested reader will find thorough instructions on how to use them in the SMC tutorial.

### Database management
`RunsFactory` is a manager dedicated to SMC runs format. All runs properties as defined from SMC argparser are automatically found.Aside from its constructor which only takes a _root_ as a single argument (path to database file), it has three main methods:
- `load(`_root_`)`: is called at the constructor to directly load database located at _root_ (if exists)
- `insert(`_rootdir_, _idx_=__-1__`)`: after the constructor, insert a SMC run located at _rootdir_ at index _idx_ in the database. By default just appends at the end of database.
- `show()`: pretty-prints the contents of the database

### Observables averaging
`Observables` is a minimal object whose main purpose is to compute time-dependant ensemble averages over a database. As a child of `RunsFactory`, its constructor requires an _ensemble_ (by default __'all'__, if not a list of indices, a slice, etc.) on top of the _root_. It has only two methods:
- `set_ensemble(`_idx_=__'all'__`)`: called at the constructor to set the ensemble. Can be re-called after initialization.
- `compute_average(`_obs_, _log_=__True__`)`:  compute the time-dependent average of _obs_ ($\in$ ["U", "MSD", "Cb", "Fs"]) over pre-defined ensemble. Of course, _obs_ must have been calculated in the SMC run.

### Plotting toolbox
`PlotToolBox` provides functions to obtain both static and dynamic configurations. Because a child of `RunsFactory`, it understands SMC runs format so that relevant configuration files are automatically read. Despite holding several methods, only two are of importance for the user:

- `render_stuff`: renders either a snapshot of particles and/or the displacement vectors from $t=0$
- `movie`: same as `render_stuff` but now dynamically

Because both of them (plus the constructor) encompass multiple user-defined specs (canvas, plot types, plot params, etc.), the user can get informations on the methods input-arguments using `help(PlotToolBox.func)` - with `func` being a specific method (`__init__` for the constructor). <br>
The tutorial remains your best assistant to understand the toolbox usage.