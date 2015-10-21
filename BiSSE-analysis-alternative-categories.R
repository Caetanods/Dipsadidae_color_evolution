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

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.

download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
load("./results_bisse.RData")

###########################

###########################
## Check for convergence using the coda package:
library(coda)

## Prepare arguments:
path <- "./mcmc_coda_plos"
burn <- 0.5
run.one <- paste("onerate_", seq(1:length(mcmc.onerate)), sep = "")
run.two <- paste("tworate_", seq(1:length(mcmc.tworate)), sep = "")

## Make trace and profile plots:
lapply(1:length(mcmc.onerate), FUN = function(x) to.mcmc.plot(mcmc.onerate[[x]][[1]], run.one[x], dir = path, burn)
lapply(1:length(mcmc.tworate), FUN = function(x) to.mcmc.plot(mcmc.tworate[[x]][[1]], run.two[x], dir = path, burn)

############################
## Run diagnostic for convergence.

hei.one <- to.heidel.diag(mcmc.onerate)
hei.two <- to.heidel.diag(mcmc.tworate)

## Summary stats of the mcmc run:
## hei.one[[1]][[1]]
## hei.two[[1]][[1]]

## Check the results of the Heidelberger and Welch's convergence diagnostic for each run:
## hei.one[[1]][[2]]
## hei.two[[1]][[2]]

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

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.

download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
load("./results_bisse.RData")

###########################

###########################
## Check for convergence using the coda package:
library(coda)

## Prepare arguments:
path <- "./mcmc_coda_plos"
burn <- 0.5
run.one <- paste("onerate_", seq(1:length(mcmc.onerate)), sep = "")
run.two <- paste("tworate_", seq(1:length(mcmc.tworate)), sep = "")

## Make trace and profile plots:
lapply(1:length(mcmc.onerate), FUN = function(x) to.mcmc.plot(mcmc.onerate[[x]][[1]], run.one[x], dir = path, burn)
lapply(1:length(mcmc.tworate), FUN = function(x) to.mcmc.plot(mcmc.tworate[[x]][[1]], run.two[x], dir = path, burn)

############################
## Run diagnostic for convergence.

hei.one <- to.heidel.diag(mcmc.onerate)
hei.two <- to.heidel.diag(mcmc.tworate)

## Summary stats of the mcmc run:
## hei.one[[1]][[1]]
## hei.two[[1]][[1]]

## Check the results of the Heidelberger and Welch's convergence diagnostic for each run:
## hei.one[[1]][[2]]
## hei.two[[1]][[2]]
