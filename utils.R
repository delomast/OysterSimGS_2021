# utility functions used in other script(s)

library(dplyr)
library(stringr)
library(optiSel)

#' create G matrix with specified base pop frequencies
#' first method of VanRaden (2008), also Endelman and Jannink (2012)
#' @param g genotypes with rows inds, cols loci, coded as 0,1,2, no missing data allowed
#' @param af base generation frequencies of the non reference allele (the allele whose dosage the number corresponds to)
#'   in the same order as the columns of `genos`
createG <- function(g, af){
	if(length(af) != ncol(g)) stop("wrong af length for number of loci")
	# P matrix in VanRaden (2008)
	p <- matrix(2*(af - 0.5), nrow = nrow(g), ncol = ncol(g), byrow = TRUE)
	return(
		# subtract 1 to convert to -1,0-1 creating M,
		# then subtract P from M to create Z
		# then ZZ' and divide by P to scale G analagous to A
		tcrossprod((g - 1) - p) / (2 * sum(af *(1 - af)))
	)
}


#' function to apply full-sib family testing with an AlphaSimR object
#' pulls ID's and phenotypes for specified proportion (rounded) of individuals from each
#' full-sib family
#' 
#' @param fam an AlphaSimR "Pop-class" object containing full-sib families
#' @param propTest the proportion of each family to phenotype. One value applied
#'   to all families. The actual number of indiviudals to phenotype will be rounded.
sibTestEqual <- function(fam, propTest){
	fs <- data.frame(mother = fam@mother, father = fam@father) %>% 
		count(mother, father) %>% mutate(nPheno = round(propTest * n))
	phenos <- data.frame()
	for(i in 1:nrow(fs)){
		# select id
		toPheno <- sample(fam@id[fam@mother == fs$mother[i] & fam@father == fs$father[i]], size = fs$nPheno[i])
		temp <- cbind(toPheno, fs$mother[i], fs$father[i], as.data.frame(fam@pheno[match(toPheno, fam@id),]))
		colnames(temp) <- c("id", "mother", "father", paste0("Trait_", 1:fam@nTraits))
		phenos <- phenos %>% rbind(temp)
	}
	return(phenos)
}

#' simple heuristic to choose markers for a genetic panel
#' this is mainly for prototyping, more sophisticated and quicker (at large scales) 
#' processes should be used
#' 
#' moves a window step wise along each chromosome and picks the highest variability marker in that window.
#' window size is calculated as numLoci/genetic length.
#' equal number of SNPs per chromosome.
#' some rounding may affect final number of SNPs.
#' 
#' @param num number of loci desired in the panel
#' @param genos genotypes (output of `AlphaSimR::pullSnpGeno`)
#' @param map genetic map of loci (output of `AlphaSimR::getSnpMap`)
quickChooseLoci <- function(num, genos, map){
	He <- (colSums(genos) / (2 * nrow(genos))) %>% (function(x) 2 * x * (1-x))() # H_exp
	map <- map %>% left_join(data.frame(He = He, id = names(He)), by = "id")
	# calculate window sizes
	distr <- map %>% group_by(chr) %>% summarise(l = max(pos)) %>% 
		mutate(num = round(num / n_distinct(chr)), wSize = l / num)
	panel <- data.frame()
	for(i in 1:nrow(distr)){ # for each chr
		skipBool <- FALSE
		cands <- map %>% filter(chr == distr$chr[i])
		for(j in 1:distr$num[i]){ # for each desired SNP
			start <- (j-1) * distr$wSize[i]
			end <- start + distr$wSize[i]
			while(TRUE){
				chosen <- cands %>% filter(pos > start & pos <= end) %>% slice(which.max(He))
				if(nrow(chosen) < 1){
					# expand window if no SNPs are available
					start <- start - distr$wSize[i] / 10
					end <- end + distr$wSize[i] / 10
					if(start < 0 && end > distr$l[i]){
						warning("ran out of SNPs in chromosome ", distr$chr[i])
						skipBool <- TRUE
						break
					}
				} else {
					# add SNP to panel
					panel <- panel %>% bind_rows(chosen)
					break
				}
			}
			if(skipBool) break
		}
	}
	return(panel)
}

#' Score function for greedy algorithm
#' from Matukumalli et al. 2009 https://doi.org/10.1371/journal.pone.0005350
#' @param start start of interval
#' @param end end of interval
#' @param maf maf of markers within interval
#' @param pos position of marker within interval (order
#'   corresponds to order of `maf`)
scoreGreedy <- function(start, end, maf, pos){
	return(
		maf * (end - start - abs((2 * pos) - end - start))
	)
}


# num = snpMap %>% select(chr) %>% distinct %>% mutate(num = 1000)

#' greedy algorithm to choose markers for a genetic panel
#' from Matukumalli et al. 2009 https://doi.org/10.1371/journal.pone.0005350
#' Assumes there are enough SNPs in each chromosome to meet the requested number
#' Note that it will not choose the first or last marker on the chromosome until
#' the very end
#' 
#' 
#' @param num data.frame with chromosome name (in map) in the first column and
#'   number of loci to choose in the second column
#' @param genos genotypes (output of `AlphaSimR::pullSnpGeno`)
#' @param map genetic map of loci (output of `AlphaSimR::getSnpMap`)
greedyChooseLoci <- function(num, genos, map){
	colnames(num) <- c("chr", "num")
	# calc minor allele freq
	maf <- (colSums(genos) / (2 * nrow(genos)))
	maf[maf > 0.5] <- 1 - maf[maf > 0.5]
	map <- map %>% left_join(data.frame(maf = maf, id = names(maf)), by = "id") 

	panel <- data.frame()
	for(i in 1:nrow(num)){ # for each chr
		cands <- map %>% filter(chr == num$chr[i])
		if(nrow(cands) < num$num[i]){
			warning("Not enough SNPs in chromosome ", num$chr[i])
			break
		}
		# calculate scores at the start
		cands <- cands %>% 
			mutate(score = scoreGreedy(start = 0, end = max(.[["pos"]]), maf = maf, pos = pos),
						 lastStart = 0,
						 lastEnd = max(.[["pos"]]))
		numSNPs <- 0
		while(TRUE){ # for each desired SNP
			# choose SNP
			temp <- which.max(cands$score)
			chosen <- cands %>% slice(temp)
			panel <- panel %>% bind_rows(chosen)
			cands <- cands[-temp,]
			numSNPs <- numSNPs + 1
			if(numSNPs == num$num[i]) break

			# recalculate scores for affected SNPs
			# only those whose interval used for the last calculation are affected
			# note that there is only one interval to update b/c the new SNP can only
			# have been in one interval
			toUpdate <- cands %>% filter(lastStart < chosen$pos & lastEnd > chosen$pos) %>%
				select(lastStart, lastEnd) %>% distinct() # this is some bs to help vectorize operations for R
			tempBool <- cands$lastStart == toUpdate$lastStart & cands$pos < chosen$pos
			cands$score[tempBool] <- scoreGreedy(start = toUpdate$lastStart,
																					 end = chosen$pos,
																					 maf = cands$maf[tempBool],
																					 pos = cands$pos[tempBool])
			cands$lastEnd[tempBool] <- chosen$pos
			tempBool <- cands$lastStart == toUpdate$lastStart & cands$pos > chosen$pos
			cands$score[tempBool] <- scoreGreedy(start = chosen$pos,
																					 end = toUpdate$lastEnd,
																					 maf = cands$maf[tempBool],
																					 pos = cands$pos[tempBool])
			cands$lastStart[tempBool] <- chosen$pos
		}
	}
	return(panel %>% arrange(chr, site))
}



#' random choice of SNPs
#' 
#' @param num number of loci desired in the panel
#' @param genos genotypes (output of `AlphaSimR::pullSnpGeno`)
#' @param map genetic map of loci (output of `AlphaSimR::getSnpMap`)
randChooseLoci <- function(num, genos, map){
	varSNPs <- colnames(genos)[!colSums(genos) %in% c(0,nrow(genos) * 2)]
	map <- map[map$id %in% varSNPs,] # select only from variable SNPs
	panel <- map[sort(sample(1:nrow(map), num, replace = FALSE)),]
	return(panel)
}


#' read in a VCF to load into AlphaSimR
#' splits the "chromosome" into separate chromosomes
#' maintains phase, assumes diploidy, keeps ony biallelic loci
#' assumes no missing genotypes (these are being treated as "true" to start a
#' simulation)
#' outputs as a matrix with rows haplotypes and columns loci
#' @param vcfPath path to the input vcf file
#' @param numLines read and evaluate in chunks of this size to reduce memory usage a little
#' @return a matrix with rows as haplotypes (adjacent rows are individuals) and cols as loci
read_vcf_for_AlphaSimR <- function(vcfPath, numLines = 20000){
	
	f <- file(vcfPath, "r") # open vcf
	on.exit(close(f))
	# move to end of header
	line <- readLines(f, n = 1)
	while(length(line) > 0){
		if(substr(line, 1, 6) == "#CHROM") break
		line <- readLines(f, n = 1)
	}
	# read in individual names
	indNames <- str_split(line, "\t")[[1]][-(1:9)]
	
	allGenos <- matrix(nrow = length(indNames) * 2, ncol = 0) # storage for all genotypes
	# now read in chunks
	genos <- readLines(f, n = numLines)
	while(length(genos) > 0){
		# splits columns by tabs and genotypes by "|"
		# so starting with col 10, each col is a haplotype and each pair
		# of columns is an individual (i.e. col 10 and 11, col 12 and 13, ...)
		genos <- str_split(genos, "\t|\\|", simplify = TRUE)
		
		locusNames <- paste0(genos[,1], "_", genos[,2]) # use CHR_position as locus names

		# now we select only the genotypes and convert to numeric, and transpose
		# to get haplotypes as rows and loci as columns
		genos <- t(matrix(as.numeric(genos[,10:ncol(genos)]), ncol = (ncol(genos) - 9)))
		colnames(genos) <- locusNames
		
		# remove any with more than 2 alleles 
		# (note: assumes no NA values, also removes any with non-standard coding)
		if(any(!(genos %in% c(0,1)))){
			# only spend time looking locus by locus if there is a hit
			toRemove <- apply(genos, 2, function(x) any(!(x %in% c(0,1))))
			genos <- genos[,!toRemove][1:10,1:10]
		}
		
		allGenos <- cbind(allGenos, genos)
		
		genos <- readLines(f, n = numLines)	
	}

	# add ind names as rows
	temp_rn <- rep(NA, nrow(allGenos))
	for(i in 1:length(indNames)){
		temp_rn[(i*2) - 1] <- paste0(indNames[i], "_hap1")
		temp_rn[i*2] <- paste0(indNames[i], "_hap2")
	}
	rownames(allGenos) <- temp_rn
	
	return(allGenos)
}

#' using results of `opticont` from `optiSel` package
#' calculate number of matings for each individual with proper rounding
#' distributes rounding error randomly according to weights
#' but with no individual having more than 1 randomly allocated
#' to it
#' @param ocsParent $parent componenet of optisolve result
#' @param N number of matings to perform
calcNumMatings <- function(ocsParent, N){
	# males
	males <- ocsParent %>% filter(Sex == "male") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(males$n)
	if(diff != 0){
		temp <- sample(1:nrow(males), size = abs(diff), prob = males$n, replace = FALSE)
		if(diff > 0) males$n[temp] <- males$n[temp] + 1
		if(diff < 0) males$n[temp] <- males$n[temp] - 1
	}
	# females
	females <- ocsParent %>% filter(Sex == "female") %>% 
		mutate(n = round(2 * N * oc))
	diff <- N - sum(females$n)
	if(diff != 0){
		temp <- sample(1:nrow(females), size = abs(diff), prob = females$n, replace = FALSE)
		if(diff > 0) females$n[temp] <- females$n[temp] + 1
		if(diff < 0) females$n[temp] <- females$n[temp] - 1
	}
	
	return(rbind(males, females))
}

#' Calculate the maximum mean kinship corresponding to 
#' a given effective population size and starting mean
#' kinship
#' @param kBar mean kinship at time 0
#' @param Ne desired effective population size
#' @param t0 start time
#' @param t end time
#' @param L generation length (time)
ubKin <- function(kBar, Ne, t0 = 0, t = 1, L = 1){
	return(1 - ((1 - kBar)*((1 - (1 / (2*Ne)))^((t - t0)/L))))
}

#' choose matings with OCS followed by inbreeding 
#' minimization
#' @param ocsData a data frame with each row corresponding to
#'   a selection candidate. Columns are Indiv, Sex (male and female)
#'   and gebv
#' @param Gmat genomic relationship matrix (the function internally converts
#'   this to the coancestry matrix)
#' @param N the number of matings to be performed
#' @param Ne the desired minimum effective population size
runOCS <- function(ocsData, Gmat, N, Ne = 50){
	# convert to coancestry/kinship matrix
	# and making sure order/presence of individuals is correct
	Gmat <- Gmat[ocsData$Indiv,ocsData$Indiv] / 2
	# data processing
	ocsCandes <- candes(phen = ocsData, N = N * 2, kin = Gmat)
	# optimum contributions
	ocsContrib <- opticont(method = "max.gebv", cand = ocsCandes, 
												 con = list(ub.kin = ubKin(kBar = ocsCandes$mean$kin, Ne = Ne)), 
												 trace=FALSE)
	# calculate number of matings per individual from contribution proportions
	ocsMatings <- calcNumMatings(ocsParent = ocsContrib$parent, N = N)
	# assign crosses (limit each pair to one cross)
	# to minimize inbreeding of each family
	crosses <- matings(ocsMatings, Kin=Amat[ocsMatings$Indiv,ocsMatings$Indiv], ub.n = 1)
	return(crosses)
}
