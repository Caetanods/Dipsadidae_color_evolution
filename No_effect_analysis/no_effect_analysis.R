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
res.full <- lapply(1:100, function(x) sim.character(phy[[x]], as.numeric(mn.full[,x]), x0=1, model = "mk2") )
sim.full.null <- t(do.call(rbind, res.full))
res.cons <- lapply(1:100, function(x) sim.character(phy[[x]], as.numeric(mn.cons[,x]), x0=1, model = "mk2") )
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

save(full.unres, full.st.bisse, cons.st.bisse, cons.unres, file = "no_effect_sim_data.RData")

## Run a MLE BiSSE analysis with all the 100 trees:
## Both with the full and cons simulated datasets.

library(multicore)
index <- 1:100

mle.cons.no.effect <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[x]], st = cons.st.bisse[[x]]
          , unres = cons.unres[[x]], constrain = "FALSE", flag = paste("no_effect_mle", x, sep=""))
          , mc.cores = 20)
save(mle.cons.no.effect, file = "no_effect_cons_mle_results.RData")

## mle.full.no.effect <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tr[[x]], st = full.st.bisse[[x]]
##           , unres = full.unres[[x]], constrain = "FALSE", flag = paste("no_effect_mle", x, sep=""))
##           , mc.cores = 20)
## save(mle.full.no.effect, file = "no_effect_full_mle_results.RData")
