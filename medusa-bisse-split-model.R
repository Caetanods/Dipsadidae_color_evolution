## This analysis use MEDUSA to search for nodes where net diversification rates change.
## Using the MEDUSA information to fit a split BiSSE model to check the influence of heterogeinity
##   in the paramter estimates.

library(multicore)
library(diversitree)
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
shifts <- lapply(1:4, FUN = function(x) medusa(tree.genus[[x]], richness = rich, criterion = "aicc", partitions = 2) )

## Prepare nodes vector for BiSSE.split
nodes <- sapply(1:4, FUN = function(x) shifts[[x]]$summary$Shift.Node[2])
nodes <- as.numeric(nodes)

## Prepare states vector for BiSSE
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Making an initial analysis only with 4 trees.
## Maybe will need more steps in order to converge.
tasks <- list(
    job1 <- function() run.bisse.split(tree = tree.genus[[1]], st = state, unres = unres, nodes[1], tun.steps = 100, chain.steps = 10000, constrain = "FALSE", flag = "phy_1"),
    job2 <- function() run.bisse.split(tree = tree.genus[[2]], st = state, unres = unres, nodes[2], tun.steps = 100, chain.steps = 10000, constrain = "FALSE", flag = "phy_2"),
    job3 <- function() run.bisse.split(tree = tree.genus[[3]], st = state, unres = unres, nodes[3], tun.steps = 100, chain.steps = 10000, constrain = "FALSE", flag = "phy_3"),
    job4 <- function() run.bisse.split(tree = tree.genus[[4]], st = state, unres = unres, nodes[4], tun.steps = 100, chain.steps = 10000, constrain = "FALSE", flag = "phy_4")
)

out <- mclapply(tasks, function(f) f(), mc.cores = 4)

save(out, file = "split_model_test.RData")
