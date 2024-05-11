#!/bin/bash

#SBATCH --output=/global/homes/h/hanyuz/AbacusHODGuide/logs/%x_%j.log
#SBATCH --error=/global/homes/h/hanyuz/AbacusHODGuide/logs/%x_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --time=30:00
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --account=desi

# Load necessary modules
module load python
conda deactivate
conda activate hod_guide
module load cray-mpich-abi

cd /global/homes/h/hanyuz/AbacusHODGuide
srun python -m mpi4py.futures scripts/run_nautilus.py --path2config configs/${JOB_NAME}.yaml
