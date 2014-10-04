## This is script is going over some trees to check if breaking the phylo into subtrees using medusa
## can improve the parameter estimates.
## Important to note that the number of parameters will increase 100% to each break in the phylogeny.
## There is a possibility that the BiSSE mcmc will not even start.

library(multicore)
library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

total <- unres[,c(1,4)]
single <- tree.genus[[1]]$tip.label[which(!tree.genus[[1]]$tip.label %in% total$tip.label)]
single <- cbind(single,1)
colnames(single) <- colnames(total) <- c("taxon","n.taxa")
rich <- rbind(total, single)

shifts <- lapply(1:5, FUN = function(x) medusa(tree.genus[[x]], richness = rich, criterion = "aicc", partitions = 2) )
nodes <- sapply(1:5, FUN = function(x) shifts[[x]]$summary$Shift.Node[2])

state <- as.numeric(st[,2])
names(state) <- st[,1]

tasks <- list(
    job1 <- function() run.bisse.split(tree = tree.genus[[1]], st = state, unres = unres, nodes[1], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"),
    job2 <- function() run.bisse.split(tree = tree.genus[[2]], st = state, unres = unres, nodes[2], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"),
    job3 <- function() run.bisse.split(tree = tree.genus[[3]], st = state, unres = unres, nodes[3], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"),
    job4 <- function() run.bisse.split(tree = tree.genus[[4]], st = state, unres = unres, nodes[4], tun.steps = 100, chain.steps = 10000, constrain = "FALSE")
)

out <- mclapply(tasks, function(f) f(), mc.cores = 4)

save(out, file = "medusa_check.RData")
