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

download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
load("./results_bisse.RData")

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

################################
## The combined posterior with half burn percentage.
## And the graphs to summarize the posterior for the parameters.

burn.comb <- list()
for(i in 1:10){
    burn.comb[[i]] <- res.cr[[i]][[1]][(dim(res.cr[[i]][[1]])[1]*50/100):dim(res.cr[[i]][[1]])[1],]
}
comb <- do.call(rbind, burn.comb)
head(comb)
dim(comb)
comb.one.rate <- comb

## save(comb.one.rate, file = "combine.one.rate.RData")
## write.table(comb, file = "Comb.posterior.bisse.cr.txt", sep = ",")
