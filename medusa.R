## This analysis use MEDUSA to search for nodes where net diversification rates change.

library(multicore)
library(geiger)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Make richness table for MEDUSA.
total <- unres[,c(1,4)]
single <- tree.genus[[1]]$tip.label[which(!tree.genus[[1]]$tip.label %in% total$tip.label)]
single <- cbind(single,1)
colnames(single) <- colnames(total) <- c("taxon","n.taxa")
rich <- rbind(total, single)

## Constraining the MEDUSA model to find at max two partitions in the phylogenies.
shifts <- mclapply(1:100, FUN = function(x) medusa(tree.genus[[x]], richness = rich, criterion = "aicc")
                 , mc.cores = 20)

save(shifts, file = "medusa.RData")
