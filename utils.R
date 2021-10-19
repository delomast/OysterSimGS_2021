# utility functions used in other script(s)

library(tibble)
library(dplyr)

#' function to apply full-sib family testing with an AlphaSimR object
#' pulls ID's and phenotypes for specified proportion (rounded) of individuals from each
#' full-sib family
#' 
#' @param fam an AlphaSimR "Pop-class" object containing full-sib families
#' @param propTest the proportion of each family to phenotype. One value applied
#'   to all families. The actual number of indiviudals to phenotype will be rounded.
sibTestEqual <- function(fam, propTest){
	fs <- tibble(mother = fam@mother, father = fam@father) %>% 
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
	map <- map %>% left_join(tibble(He = He, id = names(He)), by = "id")
	# calculate window sizes
	distr <- map %>% group_by(chr) %>% summarise(l = max(pos)) %>% 
		mutate(num = round(num / n_distinct(chr)), wSize = l / num)
	panel <- tibble()
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
