## Here I am simulating datasets using a Mk model with no effect in the speciation. I am doing this with transition
##    rates estimated either under the full and constrained models.
## After simulation I calculate the MLE of the full and constrained models. Then make a LRT.
## I will plot the LRT of the empirical data and the distribution of the simulated datasets.
## This will be similar to the plots that Rabosky and Goldberg show in their paper.

library(diversitree)
library(geiger)
source("./functions/data-prepare.R")
source("./functions/analysis.R")

## Load trees with all the species. Topology resolved randomly in BEAST.
load("./data/resolved_BD_phylo.RData")

## Get results from FigShare.
dir.create("./mcmc_BiSSE_results")
download.file(url="http://files.figshare.com/1696849/results_100_phylo_bisse.RData"
        , destfile="./mcmc_BiSSE_results/results_100_phylo_bisse.RData")

## Load results from the BiSSE mcmc for 100 trees.
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")

## Get the mean of the transition rates for each posterior for 100 trees.
## Both results from the full and constrained BiSSE models.
mn.full <- sapply(1:100, function(x) apply(mcmc.tworate[[x]][[1]][,c(6,7)], 2, mean) )
mn.cons <- sapply(1:100, function(x) apply(mcmc.onerate[[x]][[1]][,c(3,4)], 2, mean) )

## Simulate trait data under a Mk2 model with root = 1. Root value derived from ancestral state
##     estimate using BiSSE and the empirical data and trees.
## All simulations are ONLY FOR THE FIRST TREE.
res.full <- lapply(1:100, function(x) sim.character(phy[[1]], as.numeric(mn.full[,1]), x0=1, model = "mk2") )
sim.full.null <- t(do.call(rbind, res.full))
res.cons <- lapply(1:100, function(x) sim.character(phy[[1]], as.numeric(mn.cons[,1]), x0=1, model = "mk2") )
sim.cons.null <- t(do.call(rbind, res.cons))

## Drop tips to make tree a clade tree for accounting phylogenetic uncertainty while doing the
##     BiSSE estimates.
tips <- phy[[1]]$tip.label
genera <- sapply(tips, function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = FALSE )
tip.matrix <- data.frame(genera, tips, stringsAsFactors=FALSE)
tr <- to.genus.tree(phy, tip.matrix)

## Create unresolved object to run BiSSE model.
## There is the same number of species (and states) as the original phylogeny.
## They should match without problem.

null.full.states <- data.frame(species = rownames(sim.full.null), sim.full.null, stringsAsFactors=FALSE)
gen <- sapply(null.full.states[,1], function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = FALSE )
null.full.states <- data.frame(genera = gen, null.full.states, stringsAsFactors=FALSE)
rownames(null.full.states) <- NULL

null.cons.states <- data.frame(species = rownames(sim.cons.null), sim.cons.null, stringsAsFactors=FALSE)
gen <- sapply(null.cons.states[,1], function(x) strsplit(x, split = "_")[[1]][1], USE.NAMES = FALSE )
null.cons.states <- data.frame(genera = gen, null.cons.states, stringsAsFactors=FALSE)
rownames(null.cons.states) <- NULL

## Create the list of matrices to run BiSSE.
## Custom functions in 'functions/data-prepare.R'
full.unres <- make.unres(null.full.states)
full.st.bisse <- make.bisse.states(null.full.states)
cons.unres <- make.unres(null.cons.states)
cons.st.bisse <- make.bisse.states(null.cons.states)

save(full.unres, full.st.bisse, cons.st.bisse, cons.unres, file = "no_effect_LRT_sim_data.RData")

## Run a MLE BiSSE analysis with all the 100 trees:
## Both with the full and cons simulated datasets.

library(multicore)
index <- 1:100

mle.cons.fit.cons <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[1]], st = cons.st.bisse[[x]]
          , unres = cons.unres[[x]], constrain = "TRUE", flag = paste("cons.fit.cons", x, sep=""))
          , mc.cores = 20)
save(mle.cons.fit.cons, file = "mle.cons.fit.cons.RData")
rm(mle.cons.fit.cons)

mle.cons.fit.full <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[1]], st = cons.st.bisse[[x]]
          , unres = cons.unres[[x]], constrain = "FALSE", flag = paste("cons.fit.full", x, sep=""))
          , mc.cores = 20)
save(mle.cons.fit.full, file = "mle.cons.fit.full.RData")
rm(mle.cons.fit.full)

mle.full.fit.cons <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[1]], st = full.st.bisse[[x]]
          , unres = full.unres[[x]], constrain = "TRUE", flag = paste("full.fit.cons", x, sep=""))
          , mc.cores = 20)
save(mle.full.fit.cons, file = "mle.full.fit.cons.RData")
rm(mle.full.fit.cons)

mle.full.fit.full <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[1]], st = full.st.bisse[[x]]
          , unres = full.unres[[x]], constrain = "FALSE", flag = paste("full.fit.cons", x, sep=""))
          , mc.cores = 20)
save(mle.full.fit.full, file = "mle.full.fit.full.RData")

#############################################
## Now get the MLE for the empirical dataset:
rm(list = ls())
source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

emp.fit.full <- run.bisse.mle(tree = tree.genus[[1]], st = state, unres = unres, constrain = "FALSE"
                            , flag = "Empirical - full")
save(emp.fit.full, file = "emp.fit.full.RData")

emp.fit.cons <- run.bisse.mle(tree = tree.genus[[1]], st = state, unres = unres, constrain = "TRUE"
                            , flag = "Empirical - cons")
save(emp.fit.full, file = "emp.fit.cons.RData")
