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
trans.rates.hisse <- TransMatMaker(hidden.states = TRUE) ## Full model. States 0A, 0B, 1A, 1B.
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12)) ## States 0A, 1A, and 1B.

## Create and fit the model for one tree.
pp.hisse <- hisse(phy[[1]], st, hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3), trans.rate=trans.rates.hisse)

## Save the workspace:
save.image("test.hisse.RData")
