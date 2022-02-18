#!/usr/bin/env python3
# python script to evaluate switch errors between
# phased SNP file and stacks haplotype file from v2.55 (pre 2.57)
# used to evaluate mbp data from Neil
# 
# trying a different way
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
		
		nonMatchLocus = 0
		se1 = [0] * len(masterSamp)
		se2 = [0] * len(masterSamp)
		numHet = [0] * len(masterSamp)
		# get variants from snp VCF
		# start and stop are 0-based, half-open
		for x in snpIn.fetch(contig = hapLoc.chrom, start = snpPosGenome[0] - 1, stop = snpPosGenome[-1]):
			for i in range(0, len(snpPosGenome)):
				if x.pos == snpPosGenome[i] and x.id.split(":")[0] == hapLoc.id.split(":")[0]:
					for j in range(0, len(masterSamp)):
						ind = masterSamp[j]
						# skip if missing genotype
						if hapLoc.samples[ind]["GT"][0] is None:
							continue
						a1 = x.alleles[x.samples[ind]["GT"][0]]
						a2 = x.alleles[x.samples[ind]["GT"][1]]
						if a1 == a2: # only evaluate if genotype is het
							continue
						numHet[j] += 1
						hap1 = hapLoc.alleles[hapLoc.samples[ind]["GT"][0]][i]
						hap2 = hapLoc.alleles[hapLoc.samples[ind]["GT"][1]][i]
						if a1 == hap1 and a2 == hap2:
							se2[j] += 1
						elif a1 == hap2 and a2 == hap1:
							se1[j] += 1
						else:
							print("non matching alleles")
							print(x)
							print(" ")
							print(i)
							print(" ")
							print(hapLoc)
							print(ind, " ", j)
							return
							# nonMatchLocus += 1
							# numHet -= 1
					break
		
		compareLocus = 0
		seLocus = 0
		for i in range(0, len(masterSamp)):
			if numHet[i] > 1:
				compareLocus += numHet[i]
				seLocus += min(se1[i], se2[i])

		print(hapLoc.chrom, hapLoc.pos)
		print(seLocus)
		print(compareLocus)
		# print(nonMatchLocus)
		if compareLocus > 0:
			print(seLocus / compareLocus)

if __name__ == '__main__':
	Main()
