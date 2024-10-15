# Post-Processing python package for SMC software
This sub-directory of the SMC software provides various tools to analyze finished runs.
It is grounded on a manually built database so that the user can analyze his  simulations with a global perspective.

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

### Database management

### Observables

### Plotting

## Tutorial
