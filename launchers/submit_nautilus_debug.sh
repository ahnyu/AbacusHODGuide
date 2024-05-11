#!/bin/bash

# Check if job name is provided
if [ -z "$1" ]; then
  echo "Usage: $0 jobname"
  exit 1
fi

JOB_NAME=$1

# Export JOB_NAME to be available as an environment variable
export JOB_NAME

# Submit the job
sbatch --job-name=${JOB_NAME} launch_nautilus_debug.sh

