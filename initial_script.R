# simulating GS
# the world is your oyster

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

# quick founder pop for prototyping
founderPop <- runMacs2(
		nInd = 200,
		nChr = 10,
		segSites = 5100,
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
SP$addTraitA(nQtlPerChr=100) # lot's of QTL
SP$setVarE(h2=0.3) # in the range of heritability for growth, meat yield, survival, etc
SP$setSexes("yes_rand") # at the time of breeding, all individuals will only be one sex
SP$addSnpChip(5000) # all non-QTL SNPs saved from simulation

pop <- newPop(founderPop)
snpGen <- pullSnpGeno(pop) # founder SNP genotypes

# SNP panel "design"
numLoci <- c(seq(10, 100, 10), seq(125, 500, 25), seq(750, 5000, 250), seq(6000, 50000, 2000)) # number of loci in each panel to compare

numLoci <- numLoci[1:5] # for quick testing

snpMap <- getSnpMap()

# choose markers
allPanels <- lapply(numLoci, quickChooseLoci, genos = snpGen, map = snpMap)

# initial spawning
offspring <- randCross(pop, nCrosses = 100, nProgeny = 60, balance = TRUE)

# phenotype training pop (sibs)
trainPhenos <- sibTestEqual(fam = offspring, propTest = 0.5)

# now calculate things with a bunch of different panels
res <- data.frame()
for(i in 1:length(allPanels)){
	# genotypes
	g <- pullSnpGeno(offspring)[,allPanels[[i]]$id]
	
	# calc GEBVs - no imputation
	Amat <- A.mat(g - 1) # centering and calculating G with first method of VanRaden (2008), also Endelman and Jannik (2012)
	p <- data.frame(id = rownames(Amat)) |> 
		left_join(trainPhenos |> select(id, Trait_1) |> rename(pheno = Trait_1), by = "id") # hard coded for first trait
	# predict genomic breeding values
	gebv <- kin.blup(data = p, geno = "id", pheno = "pheno", K = Amat)
	comp <- data.frame(id = offspring@id, gv = gv(offspring)) |> 
		left_join(data.frame(id = names(gebv$g), gebv = gebv$g), by = "id") |>
		left_join(p, by = "id")
	
	# quick inspection of one iteration
	# cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])
	# comp %>% ggplot() + aes(x = gebv, y = gv, color = is.na(pheno)) + geom_point() + 
	# 	geom_smooth(method = "lm")
	
	# snpBLUP <- mixed.solve(p$pheno, Z = g - 1) # implementation of SNP-BLUP (note it is much faster for panel size << number of individuals)
	
	# calc accuracy of prediction and save
	res <- res |> rbind(data.frame(panelNum = i, 
																 impute = FALSE, 
																 numLoci = ncol(g), 
																 acc = cor(comp$gv[is.na(comp$pheno)], comp$gebv[is.na(comp$pheno)])))
	
	# impute
	
	# calc imputation accuracy and save
	
	# calculate GEBVs
	
	# calc accuracy of prediction and save
	
}


# now testing imputation with deeper pedigree
# make a few generations with random crosses then run above

# impute and 



