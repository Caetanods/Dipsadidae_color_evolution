## This is the main script to reproduce the analysis.
## This will run the analysis. Make convergence checks and produce the graphs.

library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run example of analysis:
res.constrain <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis:
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:100
## mcmc.onerate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = state, unres = unres, tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = state, unres = unres, tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.

## download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
## load("./results_bisse.RData")
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")

###########################

###########################
## Check for convergence using the coda package:
## Will also try to use use bonsai here.
library(coda)

## Prepare arguments:
path <- "./mcmc_coda_plos"
burn <- 0.5
run.one <- paste("onerate_", seq(1:length(mcmc.onerate)), sep = "")
run.two <- paste("tworate_", seq(1:length(mcmc.tworate)), sep = "")

## Make trace and profile plots:
lapply(1:length(mcmc.onerate)
     , FUN = function(x) to.mcmc.plot(mcmc.onerate[[x]][[1]], run.one[x], dir = path, burn) )
lapply(1:length(mcmc.tworate)
     , FUN = function(x) to.mcmc.plot(mcmc.tworate[[x]][[1]], run.two[x], dir = path, burn) )

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

## Calculate the posterior of the residency time for both the full and constrained models:
head(mcmc.tworate[[1]][[1]])
dim(mcmc.tworate[[1]][[1]])
residency.st1 <- sapply(mcmc.tworate,
                        function(x) x[[1]][5000:10000,6] / ( x[[1]][5000:10000,6] + x[[1]][5000:10000,7] )
                        )
residency.st0 <- sapply(mcmc.tworate,
                        function(x) x[[1]][5000:10000,7] / ( x[[1]][5000:10000,6] + x[[1]][5000:10000,7] )
                        )
mn <- min( c(residency.st0, residency.st1) )
mx <- max( c(residency.st0, residency.st1) )
plot(density(residency.st1, from = mn, to = mx), xlim = c(mn,mx), col = "red"
   , xlab = "Residency time", lwd = 1.5, main = "")
lines(density(residency.st0, from = mn, to = mx), col = "grey", lwd = 1.5)
