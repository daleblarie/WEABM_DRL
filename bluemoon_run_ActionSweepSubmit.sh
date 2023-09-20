#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=32G
#SBATCH --time=30:00:00

cd ${SLURM_SUBMIT_DIR}                                    
# Executable section: echoing some Slurm data
echo "Starting sbatch script run_DRL.py at:`date`"
echo "Running host:    ${SLURMD_NODENAME}"
echo "Assigned nodes:  ${SLURM_JOB_NODELIST}"

source ~/scratch/mypy/bin/activate
echo $i
echo "WE ARE TRYING TO RUN THE SCRIPT/ IT SHOULD BE NAMED ON THE LINE BELOW"
echo run_DRLtest${1}.py
mpiexec -n 20 python3 -u run_DRL_${1}.py
