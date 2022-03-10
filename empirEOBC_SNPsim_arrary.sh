#!/bin/bash
# run simulations with scrm for
# performance of different SNP panel sizes
# written for execution on Ceres
#SBATCH --cpus-per-task=1  # ask for 1 cpu
#SBATCH --mem=45G # Maximum amount of memory this job will be given
#SBATCH --time=47:00:00 # ask that the job be allowed to run for 
#SBATCH --array=36-69%75 #specify how many jobs in the array and limit number running concurrently (e.g. 1-96%40)
#SBATCH --output=arrayScrm_%a.out # tell it where to store the output console text

echo "My SLURM_JOB_ID: " $SLURM_JOB_ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load bcftools
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
mkdir /90daydata/oyster_gs_sim/temp"$SLURM_ARRAY_TASK_ID"

# echo "Begin subsampling vcf"
# subsample VCF file
# 1% of varients, which is approx 100,000 SNPs
# bcftools view ../seq_data/east_consor.vcf | perl -nle 'BEGIN { srand($x) } if (/^#/){ print; next }; print if rand(1) < 0.01' > /90daydata/oyster_gs_sim/temp"$SLURM_ARRAY_TASK_ID"/subSamp_east_consor.vcf

# echo "End subsampling vcf"
# run simulation

# randomSeed iterationNumber TemporaryLocalStorageDirectory vcfInputPath
Rscript snp_empir_HPC.R $x $SLURM_ARRAY_TASK_ID /90daydata/oyster_gs_sim/ ../seq_data/east_consor.vcf 

# remove temp directory
rm -r /90daydata/oyster_gs_sim/temp"$SLURM_ARRAY_TASK_ID"

echo "Done with simulation"
