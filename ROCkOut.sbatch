#!/bin/bash
#SBATCH -J ROCkOut
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=10
#SBATCH --mem=10G
#SBATCH -t 12:00:00
#SBATCH -e ROCkOut.err
#SBATCH -o ROCkOut.out

cd $SLURM_SUBMIT_DIR

rockout="/PATH/TO/GIT/REPO/ROCkOut/rockout_main.py"

# Download the data
python ${rockout} download -d ${dir} -p ${pos} -n ${neg} -t 10

# Build the model
python ${rockout} build -d ${dir} -t 10

# Generate the plots
python ${rockout} refine-ni -d ${dir}

## Run Log
## sbatch --export dir=model,pos=pos.txt,neg=neg.txt /PATH/TO/GIT/REPO/ROCkOut/ROCkOut.sbatch
