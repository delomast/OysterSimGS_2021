# simulating GS
# the world is your oyster

# This script looks at one generation (founders and offspring)

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

# small for quick prototyping
qtlPerChr <- 10
neutralPerChr <- 100
nChr <- 10
nFound <- 40
nOffspringPerCross <- 60

# write AlphaPeel spec file (for one chromosome)
cat("nsnp, ", neutralPerChr, "\n", 
"inputfilepath, temp/apGeno.txt
pedigree, temp/ped.txt
outputfilepath, temp/apOut
runtype, multi
ncycles, 10
", file = "apSpec.txt", sep = "")

# quick founder pop for prototyping
founderPop <- runMacs2(
		nInd = nFound,
		nChr = nChr,
		segSites = qtlPerChr + neutralPerChr,
		Ne = 12000,
		bp = 7e+07,
		genLen = 1,
		mutRate = 2.5e-08,
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
SP$addSnpChip(neutralPerChr) # all non-QTL SNPs saved from simulation

pop <- newPop(founderPop)
snpGen <- pullSnpGeno(pop) # founder SNP genotypes

# SNP panel "design"
numLoci <- unique(c(seq(100, 500, 50), seq(750, 5000, 250), seq(6000, 50000, 2000), neutralPerChr * nChr)) # number of loci in each panel to compare
numLoci <- numLoci[numLoci <= neutralPerChr * nChr]
numLoci <- numLoci[1:5] # for quick testing

snpMap <- getSnpMap()

# choose markers
allPanels <- lapply(numLoci, quickChooseLoci, genos = snpGen, map = snpMap)

# initial spawning
offspring <- randCross(pop, nCrosses = nFound/2, nProgeny = nOffspringPerCross, balance = TRUE)

# phenotype training pop (sibs)
trainPhenos <- sibTestEqual(fam = offspring, propTest = 0.5)

# now calculate things with a bunch of different panels
gebvRes <- data.frame()
imputeRes <- data.frame()
for(i in 1:length(allPanels)){
	# genotypes
	g <- pullSnpGeno(offspring)[,allPanels[[i]]$id]
	
	# calc GEBVs - no imputation
	Amat <- A.mat(g - 1) # centering and calculating G with first method of VanRaden (2008), also Endelman and Jannik (2012)
	p <- data.frame(id = rownames(Amat)) %>% 
		left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
	# predict genomic breeding values
	gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
	comp <- data.frame(id = offspring@id, gv = gv(offspring)) %>% 
		left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") %>%
		left_join(p, by = "id")
	
	# quick inspection of one iteration
	# cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])
	# comp %>% ggplot() + aes(x = gebv, y = gv, color = is.na(pheno)) + geom_point() + 
	# 	geom_smooth(method = "lm")
	
	# calc accuracy of prediction and save
	gebvRes <- gebvRes %>% rbind(data.frame(panelNum = i, 
																 impute = FALSE, 
																 numLoci = nrow(allPanels[[i]]), 
																 acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
	
	# impute
	
	# AlphaPeel
	
	# make inputs
	offspGeno <- pullSnpGeno(offspring)
	offspGeno[,!colnames(offspGeno) %in% allPanels[[i]]$id] <- 9 # low density offspring
	apGeno <- rbind(pullSnpGeno(pop), offspGeno) # combine parent (high density) and offspring genotypes
	rm(offspGeno) # save some memory
	
	ped <- SP$pedigree[,1:2]
	write.table(ped, file = "temp/ped.txt", sep = " ", quote = FALSE, col.names = FALSE, 
							row.names = TRUE)
	# impute
	imputeDose <- data.frame(id = rownames(ped))
	# for each chromosome
	for(j in 1:nChr){
		tempCols <- colnames(apGeno)[grepl(paste0("^", j, "_"), colnames(apGeno))] # loci in chromosome j
		write.table(apGeno[,tempCols],
								file = "temp/apGeno.txt", sep = " ", quote = FALSE, col.names = FALSE, 
								row.names = TRUE)
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


	# Do we need to look at other programs? maybe first determine if imputation
	# with multi-locus peeling is poor
	# SHAPEIT2 with duoHMM?
	# Beagle?
	
	# calc imputation accuracy and save
	# although wiwll discuss that this means very little for GS (Calus et al 2014 doi:10.1017/S1751731114001803)
	trueGenos <- pullSnpGeno(offspring)
	# get only loci/individuals that were imputed
	imputeCalls <- imputeDose[imputeDose$id %in% rownames(trueGenos), !colnames(imputeDose) %in% allPanels[[i]]$id]
	rownames(imputeCalls) <- imputeCalls$id
	imputeCalls <- as.matrix(imputeCalls[,-1])
	trueGenos <- trueGenos[rownames(imputeCalls), colnames(imputeCalls)]
	
	imputeRes <- imputeRes %>% 
		rbind(data.frame(panelNum = i, 
										numLoci = nrow(allPanels[[i]]),
										# mean (across loci) correlation
										imputeAcc = mean(sapply(1:ncol(trueGenos), function(x) cor(trueGenos[,x], imputeCalls[,x])))))
	rm(imputeCalls) # save some memory
	rm(trueGenos)

	# calculate GEBVs
	# note that we are overwriting earlier variables used to calcualte GEBVs
	g <- imputeDose[imputeDose$id %in% offspring@id,] # offpsring imputed values for ALL loci (even genotyped)
	rownames(g) <- g$id
	g <- as.matrix(g[,-1])
	Amat <- A.mat(g - 1) # centering and calculating G with first method of VanRaden (2008), also Endelman and Jannik (2012)
	p <- data.frame(id = rownames(Amat)) %>% 
		left_join(trainPhenos %>% select(id, Trait_1) %>% rename(pheno = Trait_1), by = "id") # hard coded for first trait
	# predict genomic breeding values
	gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
	comp <- data.frame(id = offspring@id, gv = gv(offspring)) %>% 
		left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") %>%
		left_join(p, by = "id")
	
	# calc accuracy of prediction and save
	gebvRes <- gebvRes %>% rbind(data.frame(panelNum = i, 
																					impute = TRUE, 
																					numLoci = nrow(allPanels[[i]]), 
																					acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))

}


