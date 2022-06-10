# analyze and summarize results
# produces figures for ms
library(tidyverse)


# analyze results of snp sims

gebvResAll <- tibble()
imputeResAll <- tibble()
HeResAll <- tibble()

# choose simulation to summarise
# pattern <- "^m.+rda"
pattern <- ".+rda$"

sort(as.numeric(unique(gsub("[A-Z]+|[a-z]+|_+|numOff10_|\\.rda", "", dir("./rda/", pattern = pattern)))))
sort(unique(gsub("_[0-9]+\\.rda", "", dir("./rda/", pattern = pattern))))

for (f in dir("./rda/", pattern = pattern)){
	if(grepl("EOBC", f)) next # skipping EOBC b/c not publicly available data
	if(grepl("MH", f)) next # skipping MH b/c excluded from this study
	load(paste0("./rda/", f))
	f <- gsub("multGen_|_small", "", f) # making sim names more understandable
	i <- as.numeric(gsub("[A-Z]+|[a-z]+|_+|numOff10_|\\.rda", "", f))
	s <- gsub("_[0-9]+\\.rda", "", f)
	if(grepl("curGenOnly", f)){
		gebvRes <- gebvRes %>% filter(panelNum != max(panelNum)) # remove allGen
		He_res <- He_res %>% filter(panelNum != max(panelNum)) # remove allGen
		imputeRes <- data.frame() # just making sure there is none
	}
	gebvResAll <- gebvResAll %>%
		bind_rows(gebvRes %>% 
								mutate(iter = i, sim = s))
	if(nrow(imputeRes) > 0){
		imputeResAll <- imputeResAll %>%
			bind_rows(imputeRes %>% mutate(iter = i, sim = s))
	}
	HeResAll <- HeResAll %>%
		bind_rows(He_res %>% mutate(iter = i, sim = s))
}


# table of all imputation results for supplemental
imputeTable <- imputeResAll %>% mutate(numLoci = as.factor(numLoci), generation = as.factor(genNum)) %>%
	group_by(sim, numLoci, generation) %>% summarise(sd = sd(imputeAcc), imputeAcc = mean(imputeAcc)) %>%
	select(sim, numLoci, generation, imputeAcc, sd) %>% arrange(sim, numLoci, generation)

imputeTable %>% filter(generation == 3) %>% print(n = 100)

# table of all gebv results for supplemental
gebvTable <- gebvResAll %>% mutate(numLoci = as.factor(numLoci), generation = as.factor(genNum))  %>%
	group_by(sim, numLoci, generation, impute) %>% summarise(sd = sd(acc), acc = mean(acc)) %>%
	select(sim, numLoci, generation, impute, acc, sd) %>% arrange(sim, numLoci, generation, impute)

gebvTable %>% EFGLmh::dumpTable("All_gebvTable.txt")
imputeTable %>% EFGLmh::dumpTable("All_imputeTable.txt")

HeResAll %>% EFGLmh::dumpTable("All_he_res.txt")


####
# Graphs for ms
####

# Generation 3 imputation accuracy (y) vs panel size (x) for scrm greedy, random and 10offsp
impAcc1 <- imputeTable %>%
	filter(generation == 3, grepl("scrm", sim)) %>% 
	mutate(panelType = if_else(grepl("rand", sim), "Random", "Greedy"),
				 sim = setNames(c("50 offspring, greedy", "50 offspring, random", "10 offspring, greedy"), 
				 							 c("scrm", "scrm_rand", "scrm_numOff10"))[sim]) %>%
	ggplot() + aes(x = numLoci, y = imputeAcc, group = sim, shape = sim, color = sim) +
	geom_point(size = 3, position = position_dodge(width = .3)) + geom_line(position = position_dodge(width = .3)) + 
	ylim(0.8, 1) + theme(text = element_text(size=20), legend.position=c(.85,.2)) + 
	xlab("Low-density panel size") + ylab("Imputation accuracy") + guides(shape = guide_legend(title="Simulation"),
																														color = guide_legend(title="Simulation")) +
	geom_errorbar(aes(ymin = imputeAcc - sd, ymax = imputeAcc + sd), width = 0.6, position = position_dodge(width = .3))
ggsave("impAcc1.pdf", plot = impAcc1, width = 12, height = 8)

# Generation 3 GEBV accuracy (y) vs panel size (x) for MBP

gebvAcc1 <- gebvTable %>%
	filter(generation == 3, sim == "empir_MBP") %>% 
	mutate(impute = if_else(impute, "Imputed", "Low-density only")) %>%
	ggplot() + aes(x = numLoci, y = acc, group = impute, shape = impute, color = impute) +
	geom_point(size = 3, position = position_dodge(width = .3)) + geom_line(position = position_dodge(width = .3)) +
	ylim(0.2, 0.825) + theme(text = element_text(size=20)) + 
	xlab("Low-density panel size") + ylab("GEBV accuracy") + guides(color = guide_legend(title=""), shape = guide_legend(title="")) +
	geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd, width = c(rep(0.6, 18), 0.3)), position = position_dodge(width = .3))
ggsave("gebvAcc1.pdf", plot = gebvAcc1, width = 12, height = 8)

# Imputation accuracy over generations at different panel sizes for MBP

impAcc2 <- imputeTable %>%
	filter(sim == "empir_MBP", numLoci %in% c(100, 250, 500, 750, 2000, 10000)) %>%
	ggplot() + aes(x = generation, y = imputeAcc, group = numLoci, color = numLoci) +
	geom_point(size = 3, position = position_dodge(width = .3)) + 
	geom_line(position = position_dodge(width = .3)) + ylim(0.6, 1) + 
	theme(text = element_text(size=20),legend.position=c(.85,.2)) + 
	xlab("Generation") + ylab("Imputation accuracy") + guides(color=guide_legend(title="Low-density \npanel size")) +
	geom_errorbar(aes(ymin = imputeAcc - sd, ymax = imputeAcc + sd, width = 0.6), 
								position = position_dodge(width = .3))
ggsave("impAcc2.pdf", plot = impAcc2, width = 8, height = 8)

# GEBV accuracy over generations at different panel sizes for MBP w/ and w/o imputations
mL <- gebvTable %>% filter(sim == "empir_MBP") %>% pull(numLoci) %>% as.character() %>%
	as.numeric() %>% max
gebvAcc2 <- gebvTable %>% filter(sim == "empir_MBP") %>%
	mutate(impute = if_else(impute, "Impute", "Low-density only")) %>%
	filter(numLoci %in% c(100, 250, 500, 750, 2000, 10000, mL)) %>%
	ggplot() + aes(x = generation, y = acc, group = numLoci, color = numLoci) +
	facet_wrap(~impute) +
	geom_point(size = 3, position = position_dodge(width = .4)) + 
	geom_line(position = position_dodge(width = .4)) + 
	theme(text = element_text(size=20)) + ylim(0.2, 0.825) +
	geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd), position = position_dodge(width = .4)) +
	xlab("Generation") + ylab("GEBV accuracy") + guides(color=guide_legend(title="Loci"))
ggsave("gebvAcc2.pdf", plot = gebvAcc2, width = 12, height = 8)

# Supplemental

# graph of gebv results in generation 3
pdf("supp_gebvGraphs_gen3.pdf", width = 12, height = 8)
for(i in rev(unique(gebvTable$sim))){
	if(grepl("MH", i)) next
	temp <- gebvTable %>%
		filter(generation == 3, sim == i) %>% 
		mutate(impute = if_else(impute, "Imputed", "Low-density only"))
	gebvG_gen3 <- ggplot(temp) + aes(x = numLoci, y = acc, group = impute, shape = impute, color = impute) +
		geom_point(size = 3, position = position_dodge(width = .3)) + geom_line(position = position_dodge(width = .3)) +
		ylim(0, 1) + theme(text = element_text(size=20)) + 
		xlab("Low-density panel size") + ylab("GEBV accuracy") + guides(color = guide_legend(title=""), shape = guide_legend(title="")) +
		geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd, width = c(rep(0.6, nrow(temp) - 1), 0.3)), position = position_dodge(width = .3)) +
		ggtitle(i)
	print(gebvG_gen3)
}
dev.off()

# graph of imputation across generations
pdf("supp_imputeGraphs.pdf", width = 8, height = 8)
for(i in rev(unique(imputeTable$sim))){
	if(grepl("MH", i)) next
	print(
		imputeTable %>%
			filter(sim == i, numLoci %in% c(100, 250, 500, 750, 2000, 10000)) %>%
			ggplot() + aes(x = generation, y = imputeAcc, group = numLoci, color = numLoci) +
			geom_point(size = 3, position = position_dodge(width = .3)) + 
			geom_line(position = position_dodge(width = .3)) + ylim(0.25, 1) + theme(text = element_text(size=20)) + 
			xlab("Generation") + ylab("Imputation accuracy") + guides(color=guide_legend(title="Loci")) +
			ggtitle(i) + 
			geom_errorbar(aes(ymin = imputeAcc - sd, ymax = imputeAcc + sd), width = 0.6, position = position_dodge(width = .3))
	)
}
dev.off()

# graph of gebv across generations
pdf("supp_gebvGraphs.pdf", width = 8, height = 8)
for(i in rev(unique(gebvTable$sim))){
	if(grepl("MH", i)) next
	
	mL <- gebvTable %>% filter(sim == i) %>% pull(numLoci) %>% as.character() %>%
		as.numeric() %>% max
	
	print(
		gebvTable %>% filter(sim == i) %>%
			mutate(impute = if_else(impute, "Impute", "Low-density only")) %>% 
			filter(numLoci %in% c(100, 250, 500, 750, 2000, 10000, mL)) %>%
			ggplot() + aes(x = generation, y = acc, group = numLoci, color = numLoci) +
			facet_wrap(~impute) + 
			geom_point(size = 3, position = position_dodge(width = .4)) + 
			geom_line(position = position_dodge(width = .4)) + 
			theme(text = element_text(size=20)) + ylim(0.15, 0.9) +
			geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd), position = position_dodge(width = .4)) +
			xlab("Generation") + ylab("GEBV accuracy") + guides(color=guide_legend(title="Loci")) +
			ggtitle(i)
	)
}
dev.off()

# imputation and gebv accuracy
gebvTable %>% filter(!grepl("MH", sim)) %>%
	rename(gebvAcc = acc, gebvSD = sd) %>% left_join(
	imputeTable %>% rename(imputeSD = sd) %>% mutate(impute = TRUE)) %>%
	EFGLmh::dumpTable("supp_impute_gebv_table.txt")

heSummary <- HeResAll %>% filter(!grepl("MH", sim)) %>% group_by(sim, genNum, numLoci) %>%
	summarise(H_e = mean(He), H_eSD = sd(He), maf = mean(meanMaf), mafSD = sd(meanMaf))
EFGLmh::dumpTable(heSummary, "supp_he.txt")

# looking at heterozygosity in founders
heSummary %>% filter(genNum == 0) %>% group_by(sim) %>% filter(numLoci == max(numLoci))

# drop in mean He from gen 0 to gen 3
HeResAll %>% filter(!grepl("MH", sim)) %>% select(sim, genNum, numLoci, He, iter) %>%
	spread(genNum, He) %>% mutate(diff = `0` - `3`) %>% group_by(sim, numLoci) %>%
	summarise(mDiff = mean(diff), diffSD = sd(diff)) %>% EFGLmh::dumpTable("supp_HeDiff.txt")

# greedy vs rand
gVSrand <- gebvTable %>%
	filter(generation == 3, sim == "scrm" | sim == "scrm_rand", impute) %>% 
	mutate(sim = if_else(sim == "scrm", "Greedy", "Random")) %>% 
	ggplot() + aes(x = numLoci, y = acc, group = sim, shape = sim, color = sim) +
	geom_point(size = 3, position = position_dodge(width = .3)) + geom_line(position = position_dodge(width = .3)) +
	ylim(0.6, 0.8) + theme(text = element_text(size=20)) + 
	xlab("Low-density panel size") + ylab("GEBV accuracy") + guides(color = guide_legend(title=""), shape = guide_legend(title="")) +
	geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd, width = 0.3), position = position_dodge(width = .3))
ggsave("greedVSrand.pdf", plot = gVSrand, width = 12, height = 8)


###################
## current gen for GEBV estimation only
###################

gebvResAll_c <- tibble()
HeResAll_c <- tibble()

# choose simulation to summarise
pattern <- "^curGen.+rda"

sort(as.numeric(unique(gsub("[A-Z]+|[a-z]+|_+|numOff10_|\\.rda", "", dir("./rda/", pattern = pattern)))))
sort(unique(gsub("_[0-9]+\\.rda", "", dir("./rda/", pattern = pattern))))

for (f in dir("./rda/", pattern = pattern)){
	if(grepl("EOBC", f)) next # skipping EOBC b/c not publicly available data
	load(paste0("./rda/", f))
	f <- gsub("multGen_|_small", "", f) # making sim names more understandable
	i <- as.numeric(gsub("[A-Z]+|[a-z]+|_+|numOff10_|\\.rda", "", f))
	s <- gsub("_[0-9]+\\.rda", "", f)
	gebvResAll_c <- gebvResAll_c %>%
		bind_rows(gebvRes %>% 
								mutate(iter = i, sim = s))
	HeResAll_c <- HeResAll_c %>%
		bind_rows(He_res %>% mutate(iter = i, sim = s))
}

# mark which is control (all generations)
gebvResAll_c <- gebvResAll_c %>% mutate(allGen = ifelse(panelNum == max(panelNum), TRUE, FALSE))

gebvTable_c <- gebvResAll_c %>% mutate(numLoci = as.factor(numLoci), generation = as.factor(genNum))  %>%
	group_by(sim, numLoci, generation, allGen) %>% summarise(sd = sd(acc), acc = mean(acc)) %>%
	select(sim, numLoci, generation, allGen, acc, sd) %>% arrange(sim, numLoci, generation, allGen)

curGen_graph <- gebvTable %>% filter(sim == "scrm", !impute, numLoci %in% gebvTable_c$numLoci) %>% 
	mutate(allGen = TRUE) %>% 
	bind_rows(gebvTable_c %>% filter(!allGen)) %>%
	filter(numLoci !=  1000) %>% mutate(allGen = ifelse(allGen, "All generations", "Current only")) %>%
	ggplot() + aes(x = generation, y = acc, group = numLoci, color = numLoci) +
	facet_wrap(~allGen) +
	geom_point(size = 3, position = position_dodge(width = .2)) + 
	geom_line(position = position_dodge(width = .2)) + 
	theme(text = element_text(size=20)) + ylim(0.2, 0.825) +
	geom_errorbar(aes(ymin = acc - sd, ymax = acc + sd), position = position_dodge(width = .2)) +
	xlab("Generation") + ylab("GEBV accuracy") + guides(color=guide_legend(title="Loci"), guide_legend(order = 2))

ggsave("curGen_gebvAcc.pdf", plot = curGen_graph, width = 12, height = 8)

gebvTable %>% 
	bind_rows(gebvTable_c %>% filter(!allGen) %>% select(-allGen) %>% mutate(impute = FALSE)) %>%

	# imputation and gebv accuracy
gebvTable %>% filter(!grepl("MH", sim)) %>%
	bind_rows(gebvTable_c %>% filter(!allGen) %>% select(-allGen) %>% mutate(impute = FALSE)) %>%
	rename(gebvAcc = acc, gebvSD = sd) %>% left_join(
		imputeTable %>% rename(imputeSD = sd) %>% mutate(impute = TRUE)) %>%
	EFGLmh::dumpTable("supp_impute_gebv_table.txt")

heSummary <- HeResAll %>% filter(!grepl("MH", sim)) %>% group_by(sim, genNum, numLoci) %>%
	summarise(H_e = mean(He), H_eSD = sd(He), maf = mean(meanMaf), mafSD = sd(meanMaf))
EFGLmh::dumpTable(heSummary, "supp_he.txt")	
