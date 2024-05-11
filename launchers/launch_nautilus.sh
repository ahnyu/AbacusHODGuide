#!/bin/bash

#SBATCH --output=/home/hyzhang/y1hod/logs/%x_%j.log
#SBATCH --error=/home/hyzhang/y1hod/logs/%x_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --time=12:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=desi

# Load necessary modules
module load python
conda deactivate
conda activate hod_guide
module load cray-mpich-abi

cd /home/hyzhang/y1hod
srun python -m mpi4py.futures scripts/run_nautilus.py --path2config configs/${JOB_NAME}.yaml
