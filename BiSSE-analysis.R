## This is the main script to reproduce the analysis.
## This will run the analysis and make convergence checks using the coda package.
## The script is divided into the main analysis and three alternative categorizations.
## Please check the text of the manuscript for more details of the analysis.

library(diversitree)
source("./functions/analysis.R")

#############################################################################################
## MAIN ANALYSIS AND CONVERGENCE CHECK

load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run example of analysis:
res.constrain <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis:
## WARNING: It depends on 'parallel' and takes time to run.
## library(parallel)
## index <- 1:100
## mcmc.onerate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = state, unres = unres, tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = state, unres = unres, tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.
download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData"
            , destfile = "./data/results_100_phylo_bisse.RData", method = "wget")
load("./data/results_100_phylo_bisse.RData")
###########################

###########################
## Check for convergence using the coda package:
library(coda)
source("./functions/converge-and-plots.R")

## Prepare arguments:
path <- "./mcmc_coda_plos"
burn <- 0.5
run.one <- paste("onerate_", seq(1:length(mcmc.onerate)), sep = "")
run.two <- paste("tworate_", seq(1:length(mcmc.tworate)), sep = "")

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

#############################################################################################

#############################################################################################
## ALTERNATIVE CATEGORIZATION A: ANALYSIS AND CONVERGENCE CHECK.
## In this categorization the species with bright colors only in the venter are set to the cryptic
##      category.
load("./data/data_for_BiSSE-alt.RData")

stateA <- as.numeric(st.alt[[1]][,2])
names(stateA) <- st.alt[[1]][,1]

## Run example of analysis:
res.constrain.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

#############################################################################################
## ALTERNATIVE CATEGORIZATION B: ANALYSIS AND CONVERGENCE CHECK.
## In this categorization the all contrasting species defined based on juveniles are set to the
##      cryptic state. In this analysis only the adult coloration patterns are explored.

stateB <- as.numeric(st.alt[[2]][,2])
names(stateB) <- st.alt[[2]][,1]

## Run example of analysis:
res.constrain.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

#############################################################################################
## ALTERNATIVE CATEGORIZATION C: ANALYSIS AND CONVERGENCE CHECK.
## This time the category 1 are species that have color patterns that are similar to true coral snakes
##      and category 0 includes both species that are contrasting but do not look like coral snakes
##      and cryptic species that do not show contrasting patterns.

stateC <- as.numeric(st.alt[[3]][,2])
names(stateC) <- st.alt[[3]][,1]

## Run example of analysis:
res.constrain.C <- run.bisse(tree = tree.genus[[1]], st = stateC, unres = unres.alt[[3]], tun.steps = 2, chain.steps = 5, constrain = "TRUE")
res.free.C <- run.bisse(tree = tree.genus[[1]], st = stateC, unres = unres.alt[[3]], tun.steps = 2, chain.steps = 5, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.C <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateC, unres = unres.alt[[3]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.C <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateC, unres = unres.alt[[3]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)
