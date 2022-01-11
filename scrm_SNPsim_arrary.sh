#!/bin/bash
# run simulations with scrm for
# performance of different SNP panel sizes
# written for execution on Ceres
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=20G # Maximum amount of memory this job will be given
#SBATCH --time=24:00:00 # ask that the job be allowed to run for 
#SBATCH --array=1-2%2 #specify how many jobs in the array and limit number running concurrently (e.g. 1-96%40)
#SBATCH --output=arrayMacs_%a.out # tell it where to store the output console text

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load r/4.1.2

# check for random seeds
if [ ! -f randSeeds.txt ]; then
    echo "randSeeds.txt does not exist."
	exit 1
fi

# get random seed
x=$(cat randSeeds.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo "My random seed is: " $x

# make temp directory
mkdir "$TMPDIR"/temp"$SLURM_ARRAY_TASK_ID"

# run simulation
date
echo "begin scrm"
# simulate genotypes
# 1 number of haplotypes
# 2 Ne, 
# 3 mutation rate per base, 
# 4 expected number of recombinations per chromosome per generation,
# 5 output file prefix
# 6 random seed
bash ./sim_oyster_genos_scrm.sh 400 20000 0.00000008 1 "$TMPDIR"/temp"$SLURM_ARRAY_TASK_ID"/ $x
date
echo "end scrm"

# randomSeed iterationNumber TemporaryLocalStorageDirectory
Rscript multGen_macs_HPC_array.R $x $SLURM_ARRAY_TASK_ID $TMPDIR

echo "Done with simulation"
