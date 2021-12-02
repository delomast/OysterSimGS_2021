# simulating GS
# the world is your oyster

##
# This script looks at accuracy of and with imputation
#   over several generations using founder haplotypes generated
#   by msPrime
##

library(AlphaSimR)
library(tidyverse)
library(rrBLUP)

source("utils.R")

# parallel not used yet
# suggested strategy is develop script, then run iterations of whole
# script in parallel on an HPC
# library(parallel)
# library(parallelly)
# numCores <- availableCores()
# if(Sys.info()["sysname"] == "Windows") numCores <- 1


# inputs

nOffspringPerCross <- 50
nGenerations <- 3
nQTL <- 500 # total number of qtl, will be distributed proprtional to chromosome lengths
h2 <- 0.3 # heritability
# chromosome lengths should reflect what was simulated in msPrime and 
# be in the same order
chrom_lengths <- c(32650045, 65668440, 61752955, 77061148, 59691872,
									 98698416, 51258098, 57830854, 75944018, 104168038)
vcf_from_msPrime <- "sim_pop_1.vcf"



# import founder pop from MSprime

vcfHaplos <- read_vcf_for_AlphaSimR(vcf_from_msPrime)
# subsampling to get a set of SNPs that represents 
# the largest SNP chip plus the QTLs
qtlPerChr <- round(nQTL * chrom_lengths / sum(chrom_lengths)) # sum can be slightly different that nQTL b/c of rounding error
vcfHaplos <- vcfHaplos[,sort(sample(1:ncol(vcfHaplos), 10000 + sum(qtlPerChr), replace =  FALSE))] 

# split haplotypes into chromosomes and
# create genetic map, each 1 Morgan long
#  just making it proportional to position

pos <- as.numeric(gsub("^1_", "", colnames(vcfHaplos)))
vcfHaplo_list <- list()
genMap <- list()
brPoint <- 0
for (x in chrom_lengths){
	brPoint <- brPoint + x # cutoff point for the current chromosome
	tempBool <- pos <= brPoint
	vcfHaplo_list[[length(vcfHaplo_list) + 1]] <- vcfHaplos[,tempBool]
	genMap[[length(genMap) + 1]] <- (pos[tempBool] - (brPoint - x)) / x
	vcfHaplos <- vcfHaplos[,!tempBool] # save memory
	pos <- pos[!tempBool]
}
if(ncol(vcfHaplos) > 0) warning("potential error with chromosome lengths, leftover loci")
rm(vcfHaplos, pos, tempBool)



# set up AlphaSimR objects

founderPop <- newMapPop(genMap=genMap, haplotypes=vcfHaplo_list)
SP <- SimParam$new(founderPop)
SP$setTrackPed(isTrackPed = TRUE) # have AlphaSimR maintain pedigree records

SP$addTraitA(nQtlPerChr = qtlPerChr)
SP$setVarE(h2 = h2) # in the range of heritability for growth, meat yield, survival, etc
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex
SP$addSnpChip(nSnpPerChr = sapply(vcfHaplo_list, ncol) - qtlPerChr) # all non-QTL SNPs saved from simulation



# start simulation 

pop <- list()
pop[[1]] <- newPop(founderPop)
snpGen <- pullSnpGeno(pop[[1]]) # founder SNP genotypes

# SNP panel "design"
numLoci <- unique(c(seq(100, 500, 100), 750, 1000, 2000, 5000, 10000, 15000, 30000, 60000, ncol(snpGen))) # number of loci in each panel to compare
numLoci <- numLoci[numLoci <= sum(sapply(vcfHaplo_list, ncol), -qtlPerChr)]

snpMap <- getSnpMap()

# choose markers
# allPanels <- lapply(numLoci, quickChooseLoci, genos = snpGen, map = snpMap)
allPanels <- lapply(numLoci, randChooseLoci, genos = snpGen, map = snpMap)
# allPanels <- allPanels[1:3] # quick tests

# initial spawning
pop[[2]] <- randCross(pop[[1]], nCrosses = floor(nInd(pop[[1]])/2), nProgeny = nOffspringPerCross, balance = TRUE)

trainPhenos <- data.frame()
gebvRes <- data.frame()
imputeRes <- data.frame()
for(gen in 1:nGenerations){
	# phenotype training pop (sibs) of current generation adn add to phenotype data set
	trainPhenos <- rbind(trainPhenos, sibTestEqual(fam = pop[[gen + 1]], propTest = 0.6)) # phenotype 30, select from 20
	
	for(i in 1:length(allPanels)){
		# get all genotypes at smaller panel
		g <- pullSnpGeno(pop[[1]])[,allPanels[[i]]$id]
		for(j in 2:length(pop)) g <- rbind(g, pullSnpGeno(pop[[j]])[,allPanels[[i]]$id])

		# calc GEBVs - no imputation
		Amat <- A.mat(g - 1) # centering and calculating G with first method of VanRaden (2008), also Endelman and Jannik (2012)
		p <- data.frame(id = rownames(Amat)) %>% 
			left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
		# predict genomic breeding values
		gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
		# NOTE: only using _current_ generation to calculate accuracy of gebvs
		comp <- data.frame(id = pop[[gen + 1]]@id, gv = gv(pop[[gen + 1]])) %>% 
			left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") %>%
			left_join(p, by = "id")
		# calc accuracy of prediction and save
		gebvRes <- gebvRes %>% rbind(data.frame(genNum = gen,
																						panelNum = i, 
																						impute = FALSE, 
																						numLoci = nrow(allPanels[[i]]), 
																						acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
		
		# impute
		
		# AlphaPeel
		
		# make inputs
		ped <- SP$pedigree[,1:2]
		write.table(ped, file = "temp/ped.txt", sep = " ", quote = FALSE, col.names = FALSE, 
								row.names = TRUE)
		
		# get all genotypes
		g <- pullSnpGeno(pop[[1]])
		for(j in 2:length(pop)) g <- rbind(g, pullSnpGeno(pop[[j]]))
		trueGenos <- g[as.character(pop[[gen + 1]]@id),]
		# get id's for inds with high-density genotypes (parents)
		highDensInds <- unique(c(ped[,1], ped[,2]))
		highDensInds <- highDensInds[highDensInds > 0]
		g[!rownames(g) %in% highDensInds,!colnames(g) %in% allPanels[[i]]$id] <- 9 # low density offspring
		
		# impute
		imputeDose <- data.frame(id = rownames(ped))
		# for each chromosome
		for(j in 1:length(chrom_lengths)){
			tempCols <- colnames(g)[grepl(paste0("^", j, "_"), colnames(g))] # loci in chromosome j
			write.table(g[,tempCols],
									file = "temp/apGeno.txt", sep = " ", quote = FALSE, col.names = FALSE, 
									row.names = TRUE)
			
			# write AlphaPeel spec file (for one chromosome)
			cat("nsnp, ", length(tempCols), "\n", 
					"inputfilepath, temp/apGeno.txt
pedigree, temp/ped.txt
outputfilepath, temp/apOut
runtype, multi
ncycles, 10
", file = "apSpec.txt", sep = "")
			
			# run AlphaPeel
			if(Sys.info()["sysname"] == "Windows"){
				system2("AlphaPeel/AlphaPeel_windows.exe", args = "apSpec.txt")
			} else {
				stop("OS not set up")
			}
			# load results
			tempImputeDose <- read.table("temp/apOut.dosages")
			colnames(tempImputeDose) <- c("id", tempCols)
			imputeDose <- imputeDose %>% 
				left_join(tempImputeDose %>% mutate(id = as.character(id)), by = "id")
		}
		
		# calc imputation accuracy and save
		# although we will discuss that this means very little for GS (Calus et al 2014 doi:10.1017/S1751731114001803)
		# get only loci/individuals in the _current_ generation that were imputed
		imputeCalls <- imputeDose[imputeDose$id %in% rownames(trueGenos), !colnames(imputeDose) %in% allPanels[[i]]$id]
		rownames(imputeCalls) <- imputeCalls$id
		imputeCalls <- as.matrix(imputeCalls[,-1])
		trueGenos <- trueGenos[rownames(imputeCalls), colnames(imputeCalls)]
		
		# filter out nonvariable loci
		# loci can become fixed during the simulation
		# and occasionally imputation will return the same genotype for all individuals
		# correlation is undefined (0/0) when either variable is constant
		temp <- (apply(trueGenos, 2, n_distinct) > 1) & (apply(imputeCalls, 2, n_distinct) > 1)
		trueGenos <- trueGenos[,temp ]
		imputeCalls <- imputeCalls[,temp]
		rm(temp)
		
		imputeRes <- imputeRes %>% 
			rbind(data.frame(genNum = gen,
											 panelNum = i, 
											 numLoci = nrow(allPanels[[i]]),
											 # mean (across loci) correlation
											 imputeAcc = mean(sapply(1:ncol(trueGenos), function(x) cor(trueGenos[,x], imputeCalls[,x]))),
											 numLociVar = ncol(trueGenos)))
		rm(imputeCalls) # save some memory
		rm(trueGenos)
		
		# calculate GEBVs
		# note that we are overwriting earlier variables used to calcualte GEBVs
		g <- imputeDose # imputed values for ALL loci (even genotyped)
		rownames(g) <- g$id
		g <- as.matrix(g[,-1])
		Amat <- A.mat(g - 1) # centering and calculating G with first method of VanRaden (2008), also Endelman and Jannik (2012)
		p <- data.frame(id = rownames(Amat)) %>% 
			left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
		# predict genomic breeding values
		gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
		comp <- data.frame(id = pop[[gen + 1]]@id, gv = gv(pop[[gen + 1]])) %>% 
			left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") %>%
			left_join(p, by = "id")
		
		# calc accuracy of prediction and save
		gebvRes <- gebvRes %>% rbind(data.frame(genNum = gen,
																						panelNum = i, 
																						impute = TRUE, 
																						numLoci = nrow(allPanels[[i]]), 
																						acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
	}
	
	# make next generation
	if(gen < nGenerations){
		# right now this is just within family selection
		# probably fine for the purpose of this simulation
		# but might be good to replace with something more realistic 
		brood <- comp %>% left_join(data.frame(id = pop[[gen+1]]@id,
																	sex = pop[[gen+1]]@sex,
															 mother = pop[[gen+1]]@mother,
															 father = pop[[gen+1]]@father), by = "id") %>% 
			filter(is.na(pheno)) %>% # phenotyped are not available for selection
			group_by(sex, mother, father) %>% filter(gebv == max(gebv)) %>%
			ungroup() # NOTE: b/c of use of groups, this yields a tibble
		sires <- brood %>% filter(sex == "M")
		dams <- brood %>% filter(sex == "F")
		# select crosses avoiding full sibs
		crossPlan <- matrix(nrow = nrow(dams), ncol = 2)
		crossPlan[,1] <- dams$id
		for(i in 1:nrow(dams)){
			crossPlan[i,2] <- sample(sires$id[sires$mother != dams$mother[i] & sires$father != dams$father[i]], 1)
			sires <- sires %>% filter(id != crossPlan[i,2])
		}
		
		pop[[gen + 2]] <- makeCross(pop[[gen + 1]], crossPlan = crossPlan, nProgeny = nOffspringPerCross)
	}
}

save.image("endLoopSave_multGen_msPrime.rda")

# gebvRes
# imputeRes
# 
# imputeGraph <- imputeRes %>% mutate(lociPerChrom = as.factor(numLoci / nChr), generation = as.factor(genNum)) %>%
# 	ggplot() + aes(x = generation, y = imputeAcc, group = lociPerChrom, shape = lociPerChrom, color = lociPerChrom) +
# 	geom_point(size = 3) + geom_line() + ylim(0.5, 1) + theme(text = element_text(size=20)) + ggtitle("Imputation accuracy")
# 
# gebvGraph <- gebvRes %>% mutate(lociPerChrom = as.factor(numLoci / nChr), generation = as.factor(genNum)) %>%
# 	ggplot() + aes(x = generation, y = acc, group = lociPerChrom, shape = lociPerChrom, color = lociPerChrom) +
# 	facet_wrap(~impute, nrow = 2) +
# 	geom_point(size = 3) + geom_line() + theme(text = element_text(size=20)) + ggtitle("GEBV accuracy")
# 
# ggsave("imputeGraph.pdf", plot = imputeGraph, width = 8, height = 5)
# ggsave("gebvGraph.pdf", plot = gebvGraph, width = 8, height = 10.5)
