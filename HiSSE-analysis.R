## Analysis using hidden states method.
## Method described here: Beaulieu, J. M., and B. C. O’Meara. 2015. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. bioRxiv 016386.

## This method will estimate a BiSSE model in which certain traits can have hidden traits associated with them.
## With this we can ask whether the contrasting state which has higher rates of diversification is uniformily
##      diversifing faster or if only certain clades are pulling the result. After this analysis when can
##      estimate where in the phylogenetic tree was the clades that have higher rates of diversification and
##      check their coloration patterns to elect alternative hypothesis for further testing.

## HiSSE optimizes the net turnover = lambda + mu and the extinction fraction = mu / lambda.
## The 'hisse' function set net turnover to be equal to both traits per default.

library(hisse)
library(geiger)
source("./functions/data-prepare.R")

####################################
## Load the data:
## Load trees with all the species. Topology resolved randomly in BEAST.
load("./data/resolved_BD_phylo.RData")
## Load data to run analysis.
data <- read.csv("./data/coloration_data.csv", as.is = TRUE)[,-c(1,3,5)]
states <- data$Category

####################################
## Prepare data in correct format.
## Setting categories
## Aposematic category:
states[which(states == "CON")] <- 1
states[which(states == "VEN")] <- 1
states[which(states == "CON+VEN")] <- 1

## Cryptic category:
states[which(states == "CR")] <- 0
states[which(states == "VIP")] <- 0
states[which(states == "CR+VIP")] <- 0

## All other as cryptic:
states[which(!states == "1")] <- 0

## Need to drop the NAs. Will assume that those three species are state 0 (cryptic):
states[is.na(states)] <- 0

## Make states names as sequence of "Genus_sp1", "Genus_sp2", so on.
names(states) <- number.species(data$Genus)
st <- data.frame(tips = names(states), states = states)

## Make tip.label vector as sequence of "Genus_sp1", "Genus_sp2", so on.
## For every tree in the multiphylo.
for(i in 1:length(phy)){
    tips <- sapply(strsplit(phy[[i]]$tip.label, split = "_"), function(x) x[1])
    phy[[i]]$tip.label <- number.species(tips)
}

####################################
## Prepare the HiSSE model:

## Create the hidden states matrix:
## Each is a type of model considering that we have or not hidden states:
trans.bisse <- TransMatMaker(hidden.states = FALSE) ## Equivalent of the BiSSE model.
trans.full <- TransMatMaker(hidden.states = TRUE) ## Full model. States 0A, 0B, 1A, 1B.
trans.0A.1A.1B <- ParDrop(trans.full, c(2,5,7,8,9,12)) ## States 0A, 1A, and 1B.
trans.0A.0B.1A <- ParDrop(trans.full, c(3,6,9,10,11,12)) ## States 0A, 0B, and 1A.

## Create and fit the model for one tree and with random data.

## Create random data:
## The observed probability of state 1:
prob <- 213 / (213 + 381)
## Sample the state from a binomial distribution:
ss.st <- sample(x = c(0,1), size = dim(st)[1], replace = TRUE, prob = c(1-prob,prob) )
rd.st <- st
rd.st$states <- ss.st

## Run a series of models:
tasks <- list(
    # The constrained BiSSE model:
    job1 <- function() hisse(phy[[1]], rd.st, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0)
                           , trans.rate=NULL),
    # The full BiSSE model:
    job2 <- function() hisse(phy[[1]], rd.st, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0)
                           , trans.rate=NULL),
    # The constrained 0A.1A.0B.1B HiSSE model:
    job3 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,1,1), eps.anc=c(1,1,1,1)
                           , trans.rate=trans.full),
    # The constrained 0A.1A.0B.1B, but 1B is free HiSSE model:
    job4 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2)
                           , trans.rate=trans.full),
    # The constrained 0A.1A.0B.1B, but 0B is free HiSSE model:
    job5 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,2,1), eps.anc=c(1,1,2,1)
                           , trans.rate=trans.full),
    # The constrained 0A.1A.0B.1B, but 0B,1B are separated HiSSE model:
    job6 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2)
                           , trans.rate=trans.full),
    # The constrained 0A.1A.0B.1B, but 0B and 1B are free HiSSE model:
    job7 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3)
                           , trans.rate=trans.full),
    # The full 0A.1A.0B.1B HiSSE model:
    job8 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4)
                           , trans.rate=trans.full),
    # The constrained 0A.1A.1B HiSSE model:
    job9 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,0,1), eps.anc=c(1,1,0,1)
                           , trans.rate=trans.0A.1A.1B),
    # The constrained 0A.1A.1B, but 1B is free HiSSE model:
    job10 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,0,2), eps.anc=c(1,1,0,2)
                           , trans.rate=trans.0A.1A.1B),
    # The full 0A.1A.1B HiSSE model:
    job11 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3)
                           , trans.rate=trans.0A.1A.1B),
    # The constrained 0A.1A.0B HiSSE model:
    job12 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,1,0), eps.anc=c(1,1,1,0)
                           , trans.rate=trans.0A.0B.1A),
    # The constrained 0A.1A.0B, but 0B is free HiSSE model:
    job13 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,1,2,0), eps.anc=c(1,1,2,0)
                           , trans.rate=trans.0A.0B.1A),
    # The full 0A.1A.0B HiSSE model:
    job14 <- function() hisse(phy[[1]], rd.st, hidden.states=TRUE, turnover.anc=c(1,2,3,0), eps.anc=c(1,2,3,0)
                           , trans.rate=trans.0A.0B.1A)
)

out <- mclapply(tasks, function(f) f(), mc.cores = 20)

## Save the workspace:
save.image("hisse.all_models.RData")
