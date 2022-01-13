# run AlphaImpute2
# load python, activate virutal environment
# 1 output file prefix
# 2 genotype input
# 3 pedigree input
# 4 random seed
# 5 max thread for imputation
module load python_3/3.9.2
source /project/oyster_gs_sim/ai2/bin/activate

AlphaImpute2 -out $1 -genotypes $2 -pedigree $3 -seed $4 -maxthreads $5 # -ped_only
