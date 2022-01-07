#!/bin/bash
# run simulations with macs for
# performance of different SNP panel sizes
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=20G # Maximum amount of memory this job will be given
#SBATCH --time=24:00:00 # ask that the job be allowed to run for 
#SBATCH --array=1-2%2 #specify how many jobs in the array and limit number running concurrently (e.g. 1-96%40)
#SBATCH --output=SNPdisc_RNAseq_2021/arrayGentoype_%a.out # tell it where to store the output console text

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load r/

# check for random seeds
if [ ! -f randSeeds.txt ]; then
    echo "randSeeds.txt does not exist."
	exit 1
fi

# get random seed
x=$(cat randSeeds.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "My random seed is: " $x

# run simulation
Rscript multGen_macs_HPC_array.R $x $SLURM_ARRAY_TASK_ID

echo "Done with simulation"