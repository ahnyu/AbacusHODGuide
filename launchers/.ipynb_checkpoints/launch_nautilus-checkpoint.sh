#!/bin/bash

#SBATCH --output=/home/hyzhang/y1hod/logs/%x_%j.log
#SBATCH --error=/home/hyzhang/y1hod/logs/%x_%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --mem=240GB
#SBATCH --account=rrg-wperciva

# Load necessary modules
module load StdEnv/2020
module load python
module load scipy-stack
module load mpi4py

source /home/hyzhang/generalEnv/bin/activate

# Run the Python script using mpi4py
cd /home/hyzhang/y1hod
srun python -m mpi4py.futures scripts/run_nautilus.py --path2config configs/${JOB_NAME}.yaml
