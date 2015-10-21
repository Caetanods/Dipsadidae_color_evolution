## This script will run the alternative categorization for the colors of the Dipsadidae snakes.
## There are two alternative categorizations:
## a) Species which show ventral contrasting colors as "CR".
## b) Species which only the juvenile show contrastiong colors as "CR".

library(diversitree)

## Load the data and the functions:
source("./functions/data-prepare.R")
source("./functions/analysis.R")
load("./data/data_for_BiSSE-alt.RData")

###########################
## Run the analysis:
## First alternative category A:

stateA <- as.numeric(st[,2])
names(stateA) <- st[,1]

## Run example of analysis:
res.constrain.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

## Now category B:
stateB <- as.numeric(st[,3])
names(stateB) <- st[,1]

## Run example of analysis:
res.constrain.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)
