#!/bin/bash
#SBATCH -J ROCkOut
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=10
#SBATCH --mem=50G
#SBATCH -t 18:00:00
#SBATCH -e ROCkOut.err
#SBATCH -o ROCkOut.out

cd $SLURM_SUBMIT_DIR

# load conda env for rockout
source ~/.bashrc
source activate /storage/home/hcoda1/9/rconrad6/p-ktk3-0/apps/ROCkOut

echo ' '
echo Conda env: $CONDA_PREFIX

# list input params and working directory check
echo ' '
echo Input params: ${odir} ${pos} ${neg}
echo PWD: `pwd`
echo `realpath .`
echo ' '

# path to rockout
rocker="/storage/scratch1/9/rconrad6/ROCkOut/ROCkOut/rockout_main.py"

# Download the data
echo "Running ROCker Download"
echo "Running ROCker Download" >&2
python ${rocker} download -d ${odir} -p ${pos} -n ${neg} -t 10

# Build the model
echo "Running ROCker Build"
echo "Running ROCker Build" >&2
python ${rocker} build -d ${odir} -t 10 --use_diamond

# Generate the plots
echo "Running ROCker Refine-ni"
echo "Running ROCker Refine-ni" >&2
python ${rocker} refine-ni -d ${odir}

# Collect reads to simulate test metagenome
echo "Concatenating simulated raw reads"
echo "Concatenating simulated raw reads" >&2
cat model/positive/*/raw_reads/*_250_* > simulated_raw_reads.fasta
cat model/negative/*/raw_reads/*_250_* >> simulated_raw_reads.fasta

# Run read mapping (blast/diamond alignments)
echo "Running ROCker align"
echo "Running ROCker align" >&2
python ${rockout} align -i simulated_raw_reads.fasta -d model -f score -t 10     

# Run the filter
echo "Running ROCker filter"
echo "Running ROCker filter" >&2
python ${rockout} filter -d model -f score -t 10

# Run pplacer
echo "Running ROCker pplacer"
echo "Running ROCker pplacer" >&2
python ${rockout} place -d model -f score -t 10

## Run Log
## sbatch --export odir=model,pos=pos.txt,neg=neg.txt ../../00b_sbatch/ROCkOut.sbatch
