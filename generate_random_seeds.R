# Create file containing reproducible random seeds for R functions

cmdArgs <- commandArgs(trailingOnly=TRUE)

# set "master" seed
set.seed(7)

# sample seeds
randSeeds <- sample(1:1000000, cmdArgs[1], replace = FALSE)

# write out
write.table(randSeeds, "randSeeds.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
