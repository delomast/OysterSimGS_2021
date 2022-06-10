#!/usr/bin/env python3
# python script to evaluate switch errors between
# phased SNP file and stacks haplotype file from v2.55 (pre 2.57)
# used to evaluate mbp data from Neil
# 
# Seems that genotypes from stacks for SNPs and haplotypes
# are not consistent for individual SNPs, so this is written
# to skip and snps that don't have consistent genotypes
# 

import sys
import pysam
from pysam import VariantFile

def Main():
	# defaults
	haploFile = None # -h
	snpFile = None # -s
	# get command line inputs
	flag = 1
	while flag < len(sys.argv):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == "-h":
			flag += 1
			haploFile = sys.argv[flag]
		elif sys.argv[flag] == "-s":
			flag += 1
			snpFile = sys.argv[flag]
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	if haploFile is None or snpFile is None:
		print("Error: both -h and -s must be specified")
		return
	
	hapIn = VariantFile(haploFile)
	snpIn = VariantFile(snpFile)
	
	# create master sample list
	hapSamp = [x for x in hapIn.header.samples]
	snpSamp = [x for x in snpIn.header.samples]
	masterSamp = set(hapSamp).intersection(snpSamp)
	
	# loop through loci in haplotype VCF file from stacks
	for hapLoc in hapIn.fetch():
		# get start position - int
		# this is "pos" of locus as determined by stacks, seems to be beginning of cutsite
		# start = hapLoc.pos
		# SNP positions - tuple of strings, turn into ints
		snpPos = [int(x) for x in hapLoc.info.get("snp_columns")]
		# alleles - tuple of strings
		# alleles = temp.alleles
		
		# make sure at least two variants in haplotype
		if len(snpPos) < 2:
			continue
		
		# turning snpPos into genome positions
		# need to treat - and + stacks separately
		# b/c file is stacks 2.55 (pre 2.57), snp_columns is NOT reverse complemented in "-" stacks
		if hapLoc.id[-1] == "+":
			snpPosGenome = [hapLoc.pos + x - 1 for x in snpPos]
		elif hapLoc.id[-1] == "-":
			snpPosGenome = [hapLoc.pos - x + 1 for x in snpPos]
			snpPosGenome.reverse()
		else:
			print("Error detecting orientation")
			print(hapLoc)
			return
		
		# get variants from snp VCF
		# start and stop are 0-based, half-open
		snpVars = [] # list of variant objects
		relPos = [] # list of relative position in haplotype from haplotype record, 0 based
		for x in snpIn.fetch(contig = hapLoc.chrom, start = snpPosGenome[0] - 1, stop = snpPosGenome[-1]):
			for i in range(0, len(snpPosGenome)):
				if x.pos == snpPosGenome[i] and x.id.split(":")[0] == hapLoc.id.split(:)[0]:
					snpVars += [x]
					relPos += [i]
					break

		# make sure at least two variants in snp file
		if len(snpVars) < 2:
			continue
		
		# now for each individual
		# build up variants into haplotypes at het positions
		# compare het positions to haplotype vcf
		compareLocus = 0 # number of comparisons per locus
		seLocus = 0 # number of switch errors per locus
		nonMatchLocus = 0
		for ind in masterSamp:
			# skip if missing genotype
			if hapLoc.samples[ind]["GT"][0] is None:
				continue
			# number of switch errors for both possible orientations
			se1 = 0
			se2 = 0
			# number of het genos
			numHet = 0
			# for each snp in haplotype
			for i in range(0, len(snpVars)):
				a1 = snpVars[i].alleles[snpVars[i].samples[ind]["GT"][0]]
				a2 = snpVars[i].alleles[snpVars[i].samples[ind]["GT"][1]]
				if a1 == a2: # only evaluate if genotype is het
					continue
				numHet += 1
				hap1 = hapLoc.alleles[hapLoc.samples[ind]["GT"][0]][relPos[i]]
				hap2 = hapLoc.alleles[hapLoc.samples[ind]["GT"][1]][relPos[i]]
				if a1 == hap1 and a2 == hap2:
					se2 += 1
				elif a1 == hap2 and a2 == hap1:
					se1 += 1
				else:
					# print("non matching alleles")
					# print(snpVars[i])
					# print(" ")
					# print(relPos[i])
					# print(" ")
					# print(hapLoc)
					# return
					nonMatchLocus += 1
					numHet -= 1
			
			# must be 2+ het genotypes for phasing to be non trivial
			if numHet < 2:
				continue
			
			compareLocus += numHet
			seLocus += min(se1, se2)
		
		print(hapLoc.chrom, hapLoc.pos)
		print(seLocus)
		print(compareLocus)
		print(nonMatchLocus)
		if compareLocus > 0:
			print(seLocus / compareLocus)

if __name__ == '__main__':
	Main()
