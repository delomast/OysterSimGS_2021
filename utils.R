# utility functions used in other script(s)

library(dplyr)
library(stringr)

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
#' this is mainly for prototyping, more sophisticated processes should probably be used
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
