# simulating GS
# the world is your oyster

##
# This script looks at accuracy of and with imputation
#   over several generations
#   and uses MaCS
##
# This script was used to prototype the function in 
# multGen_macs_HPC_array.R

library(AlphaSimR)
library(tidyverse)
library(rrBLUP)
library(optiSel)

source("utils.R")

# parallel not used yet
# suggested strategy is develop script, then run iterations of whole
# script in parallel on an HPC
# library(parallel)
# library(parallelly)
# numCores <- availableCores()
# if(Sys.info()["sysname"] == "Windows") numCores <- 1

qtlPerChr <- 100
neutralPerChr <- 5000 # in the maximum SNP chip
nChr <- 10
macsPop <- 1000

nFound <- 200
nOffspringPerCross <- 50
nGenerations <- 3

# simulate founder population
founderPop <- runMacs2(
	nInd = macsPop,
	nChr = nChr,
	segSites = round(1e6 / nChr), # a discovery set that you then select SNPs from
	Ne = 12000,
	bp = 7e+07,
	genLen = 1,
	mutRate = 4e-08, # higher mutation rate to reflect presumed higher rate in oysters
	histNe = c(12000, 1e+05),
	histGen = c(100, 1e+06),
	inbred = FALSE,
	split = NULL,
	ploidy = 2,
	returnCommand = FALSE,
	nThreads = NULL
)

SP <- SimParam$new(founderPop)
SP$setTrackPed(isTrackPed = TRUE) # have AlphaSimR maintain pedigree records
SP$addTraitA(nQtlPerChr = qtlPerChr)
SP$setVarE(h2 = 0.3) # in the range of heritability for growth, meat yield, survival, etc
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex
SP$addSnpChip((round(1e6 / nChr) - qtlPerChr) * nChr) # all non-QTL SNPs saved from simulation


# random breeding to give realistic LD between chromosomes
# note that some SNPs may be fixed due to drift during this process
prePop <- newPop(founderPop)
for(i in 1:9) prePop <- randCross(prePop, nCrosses = macsPop/2, nProgeny = 2, balance = TRUE)

pop <- list()
# pull founders from simulated pop while avoiding full and half sibs
pop[[1]] <- randCross(prePop, nCrosses = nFound, nProgeny = 1, balance = TRUE)
snpGen <- pullSnpGeno(pop[[1]]) # founder SNP genotypes

# number of variable SNPs in the population at the start
nVarStart <- sum(!colSums(snpGen) %in% c(0,nrow(snpGen) * 2))

# SNP panel "design"

# number of loci in each panel to compare
numLoci <- sort(unique(c(seq(100, 500, 100), 750, 1000, 2000, 5000, 10000, 15000, 25000, neutralPerChr * nChr)))
if(any(numLoci > nVarStart)) {
	warning("not enough SNPs to choose from")
	# break current iteration in parallel here
}
# numLoci <- numLoci[c(1, length(numLoci))] # for quick testing

snpMap <- getSnpMap()

# choose markers
tempNum <- lapply(numLoci, function(x){
	dfOut <- data.frame(chr = 1:nChr, num = round(x/nChr)) # all chr equal length in this sim
	# account for rounding error
	diff <- x - sum(dfOut$num)
	if(diff > 0){
		temp <- sample(1:nChr, size = diff, replace = FALSE)
		dfOut$num[temp] <- dfOut$num[temp] + 1
	} else if (diff < 0){
		temp <- sample(1:nChr, size = -diff, replace = FALSE)
		dfOut$num[temp] <- dfOut$num[temp] - 1
	}
	return(dfOut)
})
allPanels <- lapply(tempNum, greedyChooseLoci, genos = snpGen, map = snpMap)
# allPanels <- lapply(numLoci, randChooseLoci, genos = snpGen, map = snpMap) # for quick testing

# save allele freqs in base pop for each panel for calculation of G
baseAlleleFreqs <- lapply(allPanels, function(x){
	tempgenos <- snpGen[,x$id]
	return(colSums(tempgenos) / (2*nrow(tempgenos)))
})

# initial spawning
pop[[2]] <- randCross(pop[[1]], nCrosses = nFound/2, nProgeny = nOffspringPerCross, balance = TRUE)

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
		Amat <- createG(g = g, af = baseAlleleFreqs[[i]]) # G with first method of VanRaden (2008), also Endelman and Jannink (2012)
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
		ped <- SP$pedigree[,1:2] # full pedigree
		# only get inds starting with founder pop
		allInds <- c()
		for(j in 1:length(pop)) allInds <- c(allInds, pop[[j]]@id)
		ped <- ped[allInds,]
		# pretend you don't know parents of founders
		ped[pop[[1]]@id,1:2] <- 0
		write.table(ped, file = "temp/ped.txt", sep = " ", quote = FALSE, col.names = FALSE, 
								row.names = TRUE)
		
		# get all genotypes for maximum size panel - assumed at the end of allPanels
		g <- pullSnpGeno(pop[[1]])[,allPanels[[length(allPanels)]]$id]
		for(j in 2:length(pop)) g <- rbind(g, pullSnpGeno(pop[[j]])[,allPanels[[length(allPanels)]]$id])
		trueGenos <- g[as.character(pop[[gen + 1]]@id),]
		# get id's for inds with high-density genotypes (parents)
		highDensInds <- unique(c(ped[,1], ped[,2]))
		highDensInds <- highDensInds[highDensInds > 0]
		g[!rownames(g) %in% highDensInds,!colnames(g) %in% allPanels[[i]]$id] <- 9 # low density offspring
		
		# impute
		imputeDose <- data.frame(id = rownames(ped))
		# for each chromosome
		for(j in 1:nChr){
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
", file = "temp/apSpec.txt", sep = "")
			
			# run AlphaPeel
			if(Sys.info()["sysname"] == "Windows"){
				system2("AlphaPeel/AlphaPeel_windows.exe", args = "temp/apSpec.txt")
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
		Amat <- createG(g = g, af = baseAlleleFreqs[[length(allPanels)]]) # G with first method of VanRaden (2008)
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
		# OCS with lagrangian
		selCands <- comp %>% filter(is.na(pheno)) %>% pull(id)
		ocsData <- data.frame(Indiv = pop[[gen + 1]]@id, Sex = if_else(pop[[gen + 1]]@sex == "M", "male", "female")) %>%
			left_join(data.frame(Indiv = names(gebv$g), gebv = gebv$g), by = "Indiv") %>%
			filter(Indiv %in% selCands)
		matingPlan <- runOCS(ocsData = ocsData, Gmat = Amat[ocsData$Indiv,ocsData$Indiv], 
										N = nFound / 2, Ne = 50)
		
# 		# OCS with genetic algorithm implemented in AlphaMate
# 		# Creating inputs
# 		
# 		# Genomic numerator relationship matrix
# 		write.table(Amat[selCands,selCands], "temp/Amat.txt", sep = "\t", quote = FALSE, 
# 								row.names = TRUE, col.names = FALSE)
# 		# GEBVs
# 		data.frame(id = selCands) %>% left_join(comp %>% select(id, gebv), by = "id") %>%
# 			write.table("temp/gebv.txt", sep = "\t", quote = FALSE, 
# 									row.names = FALSE, col.names = FALSE)
# 		# sex
# 		data.frame(id = selCands) %>% 
# 			left_join(data.frame(id = pop[[gen+1]]@id, sex = pop[[gen+1]]@sex), by = "id") %>%
# 			mutate(sex = ifelse(sex == "M", 1, 2)) %>%
# 			write.table("temp/sex.txt", sep = "\t", quote = FALSE, 
# 									row.names = FALSE, col.names = FALSE)
# # spec file
# 		cat("NrmMatrixFile , temp/Amat.txt
# SelCriterionFile , temp/gebv.txt
# GenderFile , temp/sex.txt
# NumberOfMatings , ", nFound/2,"
# NumberOfMaleParents , ", nFound, "
# NumberOfFemaleParents , ", nFound, "
# EqualizeMaleContributions , No
# EqualizeFemaleContributions , No
# TargetMaxCriterionPct , 99
# TargetCoancestryRate , 0.01
# OutputBasename , temp/am_
# Seed , 7
# Stop
# ", file = "temp/amSpec.txt", sep = "")
# 		
# 		# run AlphaMate
# 		if(Sys.info()["sysname"] == "Windows"){
# 			system.time(system2("alphaMateWindows/AlphaMate.exe", args = "temp/amSpec.txt"))
# 		} else {
# 			stop("OS not set up")
# 		}
# 		matingPlan <- read.table("temp/am_MatingPlanModeOptTarget2.txt", header = TRUE,
# 														 colClasses = "character")
#			matingPlan <- matingPlan[,c(2,3,1)]
		
		# create next generation
		pop[[gen + 2]] <- makeCross(pop[[gen + 1]], crossPlan = as.matrix(matingPlan[,1:2]), nProgeny = nOffspringPerCross)
	}
}
save.image("endLoopSave_multGen.rda")

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
