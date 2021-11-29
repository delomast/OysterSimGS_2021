#!/usr/bin/env python
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task 72
# default memory allocation on Ceres should be fine

import msprime
import math
import random
import multiprocessing as mp
import sys
import os

# function to simulate and output
def simulate(treeSeed, mutateSeed, outputPath, rateMap):
	print("running", outputPath, "now")
	# generate tree
	ts = msprime.sim_ancestry(
		samples = 200,
		ploidy = 2,
		population_size = 12000,
		recombination_rate = rateMap,
		model = [
			msprime.DiscreteTimeWrightFisher(duration=20),
			msprime.StandardCoalescent(),
		],
		random_seed = treeSeed,
	)
	
	# add mutations
	ts = msprime.sim_mutations(ts, rate = 6e-8, # about twice the rate for humans
							   random_seed = mutateSeed)
	
	# write output
	with open(outputPath, "w") as fileOut:
		ts.write_vcf(fileOut)

def Main():

	numSims = 1000 # number of simulation iterations
	
	# setting up to simulate multiple chromosomes within each individual
  
	genetic_length = 1 # genetic length of each chromosome in Morgans
	r_break = math.log(2) # rate for crossovers between chromosomes
	chrom_lengths = [32650045, # length of chromosomes
		65668440,
		61752955,
		77061148,
		59691872,
		98698416,
		51258098,
		57830854,
		75944018,
		104168038]

	# make array of map positions
	chrom_positions = [0]
	for x in range(0, len(chrom_lengths), 1):
		chrom_positions += [chrom_positions[x] + chrom_lengths[x]]
	map_positions = [chrom_positions[0],
		chrom_positions[1]]
	for x in range(1, len(chrom_lengths), 1):
		map_positions += [chrom_positions[x] + 1, chrom_positions[x+1]]

	# calculate crossover rates for each chromosome
	r = [-math.log(1 - (genetic_length / chrom_lengths[x])) for x in range(0, len(chrom_lengths), 1)]
	# arrange rates in an array for defining the rate map
	rates = [r[0]]
	for x in range(1, len(chrom_lengths), 1):
		rates += [r_break, r[x]]

	# define RateMap object
	rate_map = msprime.RateMap(position=map_positions, rate=rates)

	# select random seeds
	random.seed(7)
	# selecting without replacement so all are unique
	treeSeeds = random.sample(range(1, 1000000001, 1), numSims) # random seeds for tree formation
	mutationSeeds = random.sample(range(1, 1000000001, 1), numSims) # random seed for mutation generation

	# get number of available threads
	sys.path.append(os.getcwd()) # necessary to add cwd to path when using slurm (since it executes a copy)
	try:
		threads = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
	except KeyError:
		threads = mp.cpu_count()
	print("detected", threads, "threads available")

	# run simulations in parallel
	process_pool = mp.Pool(processes = min(threads, numSims)) # start process pool
	for i in range(0, numSims, 1):
		process_pool.apply_async(simulate, args=(treeSeeds[i], mutationSeeds[i], "sim_pop_" + str(i + 1) + ".vcf", rate_map))
	process_pool.close()  # close the pool to new tasks
	process_pool.join()  # join the processes
	# i = 0
	# simulate(treeSeeds[i], mutationSeeds[i], "sim_pop_" + str(i + 1) + ".vcf", rate_map)

if __name__ == '__main__':
	Main()
