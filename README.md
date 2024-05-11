# AbacusHODGuide
Useful scripts for HOD fitting or optimizing using AbacusHOD

## Build a conda enviorment (assuming work on nersc)
We need a mpi4py enabled enviorment
```
module load python
conda create -n hod_guide python=3.8
conda activate hod_guide
module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
MPICC="cc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py
```
Then install all dependencies
```
conda install numpy scipy astropy numba matplotlib ipykernel
```
at where ever you want to put your abacusutils code (install the package in editable mode just in case you need to edit it)
```
git clone https://github.com/abacusorg/abacusutils.git
cd abacusutils
pip install -e .[all]
```
Nautilus sampler or any other samplers you want to use (scripts are based on Nautilus, but should be easy to modify to support other samplers)
```
conda install -c conda-forge nautilus-sampler
```
Make valid for JupyterHub
```
python -m ipykernel install --user --name hod_guide --display-name hod_guide
```

## Prep data and configuration files as needed
Check PrepDataConfig.ipynb for more details
several things
- you can either pass covmat under `config['data_params']['tracer_combos']['LRG_LRG']['path2cov']`, and set `config['fit_config_params']['mockcov'] = True`, or you can just pass covmat by `config['fit_config_params']['joint']`, there is no difference for a single tracer fitting, but for multitracer case, the first option assume no correlation between tracer combos, and the second option allow you to pass a full joint covmat. Examples use joint.
- `config['fit_params']['LRG']['some_param_name'] = [paramindex, lower bound, higher bound, gaussian mean, gaussian std]`, following https://arxiv.org/abs/2306.06314, we use a truncated normal as a prior, if you do not want to apply any prior other than flat bound, then set `config['fit_params']['LRG']['some_params_name'] = [paramindex, lower bound, higher bound]`, then change scripts/run_nautilus.py

## Update Slurm file
Check launchers/launch_nautilus.sh, update path etc.
```
#SBATCH --output=/path/to/AbacusHODGuide/logs/%x_%j.log
#SBATCH --error=/path/to/AbacusHODGuide/logs/%x_%j.err
#... change number of nodes, wall time etc if needed
cd /path/to/AbacusHODGuide
```

## Full HOD fitting
say you want a full HOD fitting run using configuration file some_config_name.yaml
```
cd /path/to/AbacusHODGuide
sh launchers/submit_nautilus.sh some_config_name

