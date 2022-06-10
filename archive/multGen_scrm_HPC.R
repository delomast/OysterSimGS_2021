# This script looks at accuracy of and with imputation
#   over several generations
#   and uses scrm to simulate the genome
# It is meant to be called by Rscript as part of a
# slurm array job

print(Sys.time())
print("begin")

.libPaths(c(.libPaths(), "/project/oyster_gs_sim/R_packages/4.1/"))
library(AlphaSimR, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(tidyverse, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(rrBLUP, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")
library(optiSel, lib.loc="/project/oyster_gs_sim/R_packages/4.1/")

source("utils.R")

cmdArgs <- commandArgs(trailingOnly=TRUE)

# Script parameters given on the command line
#' @param randSeed random seed to set R's random number generator
#' @param iterationNumber used to set unique file output names
randSeed <- cmdArgs[1]
iterationNumber <- cmdArgs[2]
localTempDir <- cmdArgs[3]

set.seed(randSeed)

chrLen <- c(32650045, 65668440, 61752955, 77061148, 59691872, 98698416, 51258098, 57830854, 75944018, 104168038)

qtlPerChr <- round(1000 * (chrLen/ sum(chrLen))) # does result in 1000, no rounding error
neutralPerChr <- round(50000 * (chrLen/ sum(chrLen))) # in the maximum SNP chip per chromosome # results in 50000
nChr <- 10

nFound <- 200 # this will be checked against the number of loaded scrm haplotypes
nOffspringPerCross <- 50
nGenerations <- 3


# load in genotypes from scrm output and set up alphaSimR simulation
print(Sys.time())
print("begin loading scrm genotypes")

haplo_list <- list()
genMap <- list()
for(i in 1:nChr){
	haplo_list[[i]] <- read_scrm_transpose_for_AlphaSimR(paste0(localTempDir, "/", "temp", iterationNumber, "/chr", i, ".txt"),
													numLoci = qtlPerChr[i] + neutralPerChr[i], min_maf = 0.05, numLines = 20000, incr = 1/chrLen[i])
	genMap[[i]] <- as.numeric(gsub("^pos_", "", colnames(haplo_list[[i]])))
	# remove simulated haplotypes after loading in
	file.remove(paste0(localTempDir, "/", "temp", iterationNumber, "/chr", i, ".txt"))
}

print(Sys.time())
print("end loading")

if(any(sapply(haplo_list, nrow) != (nFound * 2))) stop("mismatch in number of founders and number of haplotypes from scrm")

founderPop <- newMapPop(genMap=genMap, haplotypes=haplo_list)
SP <- SimParam$new(founderPop)
SP$setTrackPed(isTrackPed = TRUE) # have AlphaSimR maintain pedigree records
SP$addTraitA(nQtlPerChr = qtlPerChr)
SP$setVarE(h2 = 0.3) # in the range of heritability for growth, meat yield, survival, etc
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex
SP$addSnpChip(nSnpPerChr = neutralPerChr) # all non-QTL SNPs saved from simulation

pop <- list()
# pull founders from simulated pop while avoiding full and half sibs
pop[[1]] <- newPop(founderPop)
snpGen <- pullSnpGeno(pop[[1]]) # founder SNP genotypes

print(Sys.time())
print("begin panel design")

# SNP panel "design"

# number of loci in each panel to compare
numLoci <- sort(unique(c(seq(100, 500, 100), 750, 1000, 2000, 5000, 10000, 15000, 25000, 50000)))

snpMap <- getSnpMap()

# choose markers
tempNum <- lapply(numLoci, function(x){
	dfOut <- data.frame(chr = 1:nChr, num = round(x * (chrLen/ sum(chrLen)))) # proportional to chr length
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
# optimization opportunity:
# by applying the greedy algorithm to all panels, you are repeating computation
# could apply to the largest panel and then rank by score within each chromosome
# and take the top X to produce smaller panels
allPanels <- lapply(tempNum, greedyChooseLoci, genos = snpGen, map = snpMap)

print(Sys.time())
print("end panel design")

# save allele freqs in base pop for each panel for calculation of G
baseAlleleFreqs <- lapply(allPanels, function(x){
	tempgenos <- snpGen[,x$id]
	return(colSums(tempgenos) / (2*nrow(tempgenos)))
})

# initial spawning
pop[[2]] <- randCross(pop[[1]], nCrosses = nFound/2, nProgeny = nOffspringPerCross, balance = TRUE)

if(!dir.exists(paste0(localTempDir, "/", "temp", iterationNumber))) dir.create(paste0(localTempDir, "/", "temp", iterationNumber))
trainPhenos <- data.frame()
gebvRes <- data.frame()
imputeRes <- data.frame()
for(gen in 1:nGenerations){
	print(Sys.time())
	print(paste("begin gen: ", gen))
	# phenotype training pop (sibs) of current generation adn add to phenotype data set
	trainPhenos <- rbind(trainPhenos, sibTestEqual(fam = pop[[gen + 1]], propTest = 0.6)) # phenotype 30, select from 20
	
	for(i in 1:length(allPanels)){
		print(Sys.time())
		print(paste("begin panel: ", i))
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
		if (i == length(allPanels)) next # no need to impute if all loci are genotyped
		
		# AlphaPeel
		print(Sys.time())
		print("begin impute")
		# make inputs
		ped <- SP$pedigree[,1:2] # full pedigree
		# only get inds starting with founder pop
		allInds <- c()
		for(j in 1:length(pop)) allInds <- c(allInds, pop[[j]]@id)
		ped <- ped[allInds,]
		# pretend you don't know parents of founders
		ped[pop[[1]]@id,1:2] <- 0
		write.table(ped, file = paste0(localTempDir, "/", "temp", iterationNumber, "/ped.txt"),
								sep = " ", quote = FALSE, col.names = FALSE, row.names = TRUE)
		
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
									file = paste0(localTempDir, "/", "temp", iterationNumber, "/apGeno.txt"), 
									sep = " ", quote = FALSE, col.names = FALSE, 
									row.names = TRUE)
			# write AlphaPeel spec file (for one chromosome)
			cat("nsnp, ", length(tempCols), "\n", 
					"inputfilepath, ", localTempDir, "/","temp", iterationNumber, "/apGeno.txt
pedigree, ", localTempDir, "/", "temp", iterationNumber, "/ped.txt
outputfilepath, ", localTempDir, "/", "temp", iterationNumber, "/apOut
runtype, multi
ncycles, 10
", file = paste0(localTempDir, "/", "temp", iterationNumber, "/apSpec.txt"), sep = "")
			
			# run AlphaPeel
			system2("AlphaPeel/AlphaPeel_linux", args = paste0(localTempDir, "/", "temp", iterationNumber, "/apSpec.txt"))
			
			# load results
			tempImputeDose <- read.table(paste0(localTempDir, "/", "temp", iterationNumber, "/apOut.dosages"))
			colnames(tempImputeDose) <- c("id", tempCols)
			imputeDose <- imputeDose %>% 
				left_join(tempImputeDose %>% mutate(id = as.character(id)), by = "id")
		}
		
		print(Sys.time())
		print("end impute")
		
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
		print(Sys.time())
		print("begin ocs")
		# OCS with lagrangian
		selCands <- comp %>% filter(is.na(pheno)) %>% pull(id)
		ocsData <- data.frame(Indiv = pop[[gen + 1]]@id, Sex = if_else(pop[[gen + 1]]@sex == "M", "male", "female")) %>%
			left_join(data.frame(Indiv = names(gebv$g), gebv = gebv$g), by = "Indiv") %>%
			filter(Indiv %in% selCands)
		matingPlan <- runOCS(ocsData = ocsData, Gmat = Amat[ocsData$Indiv,ocsData$Indiv], 
												 N = nFound / 2, Ne = 50)
		print(Sys.time())
		print("end ocs")
		# create next generation
		pop[[gen + 2]] <- makeCross(pop[[gen + 1]], crossPlan = as.matrix(matingPlan[,1:2]), nProgeny = nOffspringPerCross)
	}
}
save.image(paste0("multGen_scrm_", iterationNumber, ".rda"))
# for low memory use
# save(snpGen, imputeRes, gebvRes, file = paste0("multGen_scrm_small_", iterationNumber, ".rda"))
