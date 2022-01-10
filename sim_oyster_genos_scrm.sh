# simulate oyster genomes with scrm
# single thread, one chromosome at a time
# outputs will be separate files by chromosome with
# output file names $5chr[1-10].txt
# takes positional arguments in order:
# 1 number of haplotypes
# 2 Ne, 
# 3 mutation rate per base, 
# 4 expected number of recombinations per chromosome per generation,
# 5 output file prefix
# 6 random seed

# activate conda environment containing scrm
module load miniconda
source activate /project/oyster_gs_sim/scrm

# now loop over chromosomes
chrom_lengths=(32650045 65668440 61752955 77061148 59691872 98698416 51258098 57830854 75944018 104168038)
i=0
for len in "${chrom_lengths[@]}"; do
	# increment counter
	i=$((i+1))
	# calculate 4*Ne*u*L
	t=$(awk "BEGIN {print 4*$2*$3*$len}")
	# calculate 4*Ne*r
	R=$(awk "BEGIN {print 4*$2*$4}")
	scrm $1 1 -t $t -r $R $len -l 200r -seed $6 > "$5"chr"$i".txt
done
