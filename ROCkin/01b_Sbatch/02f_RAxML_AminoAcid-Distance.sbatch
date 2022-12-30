#PBS -N RAxML-D
#PBS -l nodes=1:ppn=8
#PBS -l mem=10gb
#PBS -l walltime=24:00:00
#PBS -q inferno     
#PBS -A GT-ktk3-CODA20
#PBS -o 00a_log/02f_RAxML-D_${name}.out
#PBS -e 00a_log/02f_RAxML-D_${name}.err

#!/bin/bash
#SBATCH -J RAxML-D
#SBATCH -q YOUR_QUEUE
#SBATCH -A YOUR_ACCOUNT
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH -t 4:00:00:00
#SBATCH -o 02f_RAxML-D.out
#SBATCH -e 02f_RAxML-D.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Load the modules for your servers RAxML or RAxML must be avialable in your PATH
module load raxml/8.2.11-mva2

# set some veriables.
pre=`pwd`/${output}
out=${name}.RAxML

if [ ! -d ${output} ]; then mkdir ${output}; fi

# Run RAxML for amino acid sequence
# Bootstrapped maximum likelihood tree
#raxmlHPC-PTHREADS-AVX -s ${input} -n ${out} -w ${pre} -m PROTGAMMAAUTO -f a -N autoMRE -p 4564821 -T 8 -x 1235

# single ML distance tree
raxmlHPC-PTHREADS-AVX -s ${input} -n ${out} -w ${pre} -m PROTGAMMAAUTO -f x -p 4564821 -T 8

## Run Log ##
## sbatch --export input=fasta.aln,output=directory,name=jobname /Path/to/GitHub/repo/01b_Sbatch/02f_RAxML_AminoAcid.sbatch

