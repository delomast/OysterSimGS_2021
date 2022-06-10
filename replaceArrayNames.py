#!/usr/bin/env python3
# replacing array names with chr/pos based on mapping in reference genome

def Main():
	newContigs = {"NC_047559.1" : 55785328,
								"NC_047560.1" : 73222313,
								"NC_047561.1" : 58319100,
								"NC_047562.1" : 53127865,
								"NC_047563.1" : 73550375,
								"NC_047564.1" : 60151564,
								"NC_047565.1" : 62107823,
								"NC_047566.1" : 58462999,
								"NC_047567.1" : 37089910,
								"NC_047568.1" : 57541580}
								
	# load in new names from SAM file
	newPos = {}
	with open("../guttierez_array_bowtie2_map.sam") as samIn:
		line = samIn.readline()
		while line:
			sep = line.rstrip().split("\t")
			if sep[2] in newContigs:
				newPos[sep[0]] = sep[2:4] # old ID name is key (note ID not CHROM), value is [new name, new position]
			line = samIn.readline()
	
	with open("../kijasOrig.vcf", "r") as vcfIn, open("../kijas_chr.vcf", "w") as vcfOut:
		line = vcfIn.readline()
		while line:
			sep = line.rstrip()
			if sep.startswith("##"):
				if not sep.startswith("##contig"):
					vcfOut.write(sep + "\n")
			elif sep.startswith("#CHROM"):
				# write new contig names and lengths last
				for k in newContigs:
					vcfOut.write("##contig=<ID=" + k + ",length=" + str(newContigs[k]) + ">" + "\n")
				vcfOut.write(sep + "\n")
				# and finish header loop
				break
			else:
				print("error")
				return 0
			line = vcfIn.readline()
		# now process genotypes
		line = vcfIn.readline()
		numGenos = 0 # number written out
		while line:
			sep = line.rstrip().split("\t")
			tempNewName = newPos.get(sep[2], None)
			# replace old with new chrom and pos IF mapped successfully to a chromosome
			if tempNewName is not None:
				sep[0:2] = tempNewName
				vcfOut.write("\t".join(sep) + "\n") # write out
				numGenos += 1
			line = vcfIn.readline()
		print(numGenos)

Main()
