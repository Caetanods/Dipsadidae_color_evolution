## Test for the magnitude of the difference between the null model and the empirical results
##     of the BiSSE estimation.
## The null model is a MK2 model with the transition rates estimated from the data in the full
##     BiSSE model. In this null models there is no effect of the trait in the speciation and
##     extinction.
## To simulate the data under the null model I am using a fully resolved tree. The species were
##     inserted in a polytomy and the polytomy were randomly resolved and branch lengths were
##     calculated under a strict clock model in BEAST.

library(diversitree)
library(geiger)

## Load trees with all the species. Topology resolved randomly in BEAST.
load("./data/resolved_BD_phylo.RData")
## Load results from the BiSSE mcmc for 100 trees.
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")
rm(mcmc.onerate)

## Get the mean of the transition rates for each posterior for 100 trees.
ll <- length(mcmc.tworate)
mn <- sapply(1:ll, function(x) apply(mcmc.tworate[[x]][[1]][,c(6,7)], 2, mean) )

## Simulate trait data under a Mk2 model with root = 1. Root value derived from ancestral state
##     estimate using BiSSE and the empirical data and trees.
res <- lapply(1:ll, function(x) sim.character(phy[[x]], as.numeric(mn[,x]), x0=1, model = "mk2") )
sim.null <- t(do.call(rbind, res))

## Drop tips to make tree a clade tree for accounting phylogenetic uncertainty while doing the
##     BiSSE estimates.
tips <- phy[[1]]$tip.label
genera <- sapply(tips, function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = FALSE )
tip.matrix <- data.frame(genera, tips, stringsAsFactors=FALSE)

source("./functions/data-prepare.R")
tr <- to.genus.tree(phy, tip.matrix)

## Create unresolved object to run BiSSE model.
## There is the same number of species (and states) as the original phylogeny.
## They should match without problem.

null.states <- data.frame(species = rownames(sim.null), sim.null, stringsAsFactors=FALSE)
gen <- sapply(null.states[,1], function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = FALSE )
null.states <- data.frame(genera = gen, null.states, stringsAsFactors=FALSE)
rownames(null.states) <- NULL

## Create the list of matrices to run BiSSE.
## Custom functions in 'functions/data-prepare.R'
unres <- make.unres(null.states)
st.bisse <- make.bisse.states(null.states)

## Load functions for analysis:
source("./functions/analysis.R")

## Run example of analysis:
## res.null <- run.bisse(tree = tr[[1]], st = st.bisse[[1]], unres = unres[[1]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run analysis with 20 first trees:
## This analysis will take one week to complete. I am going to use more trees in a second
##      batch of analyses if necessary.
## WARNING: It depends on 'multicore' and takes time to run.
## WARNING: It is set to use 20 cores.

library(multicore)
index <- 1:20
mcmc.no.effect <- mclapply(index, FUN = function(x) run.bisse(tree = tr[[x]], st = st.bisse[[x]]
          , unres = unres[[x]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"
          , flag = paste("no_effect_", x, sep="")), mc.cores = 20)

## Save the full workspace to get all simulation results.
save.image(file = "no_effect_run_image.RData")

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.

## Compare the posterior of the null model with the empirical data.
