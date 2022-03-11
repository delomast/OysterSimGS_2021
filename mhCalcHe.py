#!/usr/bin/env python3
# This script is called by the microhap simulation
# R script to 1) select random SNPs to use as QTL
# 2) calculate He for windows which can then be used
# to run the greedy algorithm
# 
# command line arguments in order
# input vcf file, random seed, number of random SNPs, output prefix, window size

# load libraries
import sys
import random

#' @param genoList list of split VCF rows with only genotype columns
#'   and genotypes are all biallelic
def calcHe(genoList):
	# build string representations of haplotypes
	alleles1 = [""] * len(genoList[0])
	alleles2 = [""] * len(genoList[0])
	for g in genoList:
		alleles1 = [alleles1[i] + g[i][0] for i in range(0, len(alleles1))]
		alleles2 = [alleles2[i] + g[i][2] for i in range(0, len(alleles2))]
	# count allele frequencies
	af = {}
	for a in alleles1 + alleles2:
		af[a] = af.get(a, 0) + 1
	# calc He
	He = 1
	for a in af:
		f = af[a] / (len(alleles1) + len(alleles2))
		He -= f**2
	
	# returning as string b/c is immediately being written to file
	return str(He)


def Main():
	# command line arguments
	# input file, random seed, number of random SNPs, output prefix, window size
	inFile = sys.argv[1]
	randSeed = int(sys.argv[2])
	numRand = int(sys.argv[3])
	outPrefix = sys.argv[4]
	wS = int(sys.argv[5])
	
	#######
	# Select random SNPs, outputs chr and pos to file
	
	# set random seed
	random.seed(randSeed)
	
	# determine number of SNPs in the file
	numSNPs = 0
	chrom = []
	pos = []
	with open(inFile, "r") as fileIn:
		line = fileIn.readline()
		while line:
			if line[0] != "#":
				line = line.rstrip().split("\t")
				numSNPs += 1
				chrom += [line[0]]
				pos += [line[1]]
			line = fileIn.readline()
	
	# choose positions with first SNP being 1
	qtl = random.sample(range(1,numSNPs + 1), numRand)
	qtl.sort() # sort list
	
	# write out
	with open(outPrefix + "qtl.txt", "w") as outFile:
		for i in qtl:
			outFile.write("\t".join([str(i), chrom[i - 1], pos[i-1]]) + "\n")
	
	#######
	# Calculate He with sliding window outputs chr, pos, and He to file
	# and chosen skipping qtl
	qtl.reverse() # reversing so can use pop
	
	# run sliding window across input VCF file
	with open(inFile, "r") as fileIn, open(outPrefix + "He.txt", "w") as outFile:
		
		# skip header
		line = fileIn.readline()
		while line[0] == "#":
			line = fileIn.readline()
		
		# set up beginning of first window
		lineNum = 1 # "line number" with first SNP being line 1
		if len(qtl) > 0:
			toSkip = qtl.pop()
		else:
			toSkip = None
		while lineNum == toSkip:
			line = fileIn.readline()
			lineNum += 1
			if len(qtl) > 0:
				toSkip = qtl.pop()
			else:
				toSkip = None
		line = line.rstrip().split("\t")
		curChr = line[0] # current chromosome window is on
		cur = [int(line[1])] # current positions in the window
		genos = [line[9:]] # genotypes for current window
		windLines = [lineNum] # "line numbers" for current window
		nextR = fileIn.readline() # next available position
		lineNum += 1
		
		# loop through all possible windows
		numSnpsPerWindow = {}
		numWindows = 0
		while nextR:
			# skip if it's a random qtl SNP
			if lineNum == toSkip:
				nextR = fileIn.readline()
				lineNum += 1
				if len(qtl) > 0:
					toSkip = qtl.pop()
				else:
					toSkip = None
				continue
			
			nextR = nextR.rstrip().split("\t")
			nextR[1] = int(nextR[1])
			# determine if the next variant is within the window
			if nextR[0] == curChr and (nextR[1] - cur[0]) < wS:
				cur += [nextR[1]] # add to window
				genos += [nextR[9:]]
				windLines += [lineNum]
			else:
				# evaluate current window
				numWindows += 1
				numSnpsPerWindow[len(cur)] = numSnpsPerWindow.get(len(cur), 0) + 1

				# now calculate expHet and write to output
				outFile.write("\t".join([curChr, ",".join([str(x) for x in cur])] + [calcHe(genos)] + [",".join([str(x) for x in windLines])]) + "\n")
			
				# advance to next window
				if curChr == nextR[0]:
					# have to include the next snp, otherwise
					# new window is a subset of the previous window
					newCur = []
					newGenos = []
					newWindLines = []
					for i in range(0, len(cur)):
						if (nextR[1] - cur[i]) < wS:
							newCur += [cur[i]]
							newGenos += [genos[i]]
							newWindLines += [windLines[i]]
					cur = newCur
					genos = newGenos
					windLines = newWindLines
					# adding new SNP
					cur += [nextR[1]]
					genos += [nextR[9:]]
					windLines += [lineNum]
				else:
					curChr = nextR[0]
					cur = [nextR[1]]
					genos = [nextR[9:]]
					windLines = [lineNum]
			nextR = fileIn.readline()
			lineNum += 1
			
		# while loop will end with the last window un-evaluated
		# evaluate last window
		numWindows += 1
		numSnpsPerWindow[len(cur)] = numSnpsPerWindow.get(len(cur), 0) + 1
		outFile.write("\t".join([curChr, ",".join([str(x) for x in cur])] + [calcHe(genos)] + [",".join([str(x) for x in windLines])]) + "\n")
		
		# print some summary information
		print("number of windows: ", numWindows)
		print("Distribution of SNPs per window")
		print("SNPs ", "NumberOfWindows")
		for k in sorted(list(numSnpsPerWindow)):
			print(k, numSnpsPerWindow[k])

	

if __name__ == '__main__':
	Main()
