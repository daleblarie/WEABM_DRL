#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=32G
#SBATCH --time=30:00:00
#SBATCH --job-name=WEABM_Experiment1 --output=%x.out

cd ${SLURM_SUBMIT_DIR}                                    
# Executable section: echoing some Slurm data
echo "Starting sbatch script run_DRLtest.py at:`date`"
echo "Running host:    ${SLURMD_NODENAME}"
echo "Assigned nodes:  ${SLURM_JOB_NODELIST}"

source ~/scratch/mypy/bin/activate

mpiexec -n 20 python3 -u run_DRL.py
