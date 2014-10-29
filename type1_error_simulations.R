## This script will do some simulations to check how the type 1 error affect my results.

## Sim1: Randomizing the data on the tips.
## For this I will create a new dataset of tip data generated with the same frequency of states of the
##     real dataset.

library(multicore)
library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE with real data:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Get total of species and frequency of each state.
## Sum of the 'unres' col plus the monotypic genera.
total.sp <- sum(unres$Nc) + length(which(!is.na(state)))
st1.freq <- ( sum(unres$n1) + length(which(state == "1")) ) / total.sp
st0.freq <- ( sum(unres$n0) + length(which(state == "0")) ) / total.sp

## Function to create the 'unres' and 'state' objects from sampled data with defined frequencies:
rd.data <- to.make.rd.state(state, unres, st0.freq, st1.freq)

## Run example of analysis under the full model.
#res.constrain <- run.bisse(tree = tree.genus[[1]], st = rd.data$random.st, unres = rd.data$random.unres
#                         , tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Sim2: Simulating a no-effect MK2 model.
## For this we need the complete resolved phylogeny at the species level.

## First use phytools to create the polytomies:
library(phytools)

## Create the diver object. Equal to the rich parameter for MEDUSA.
diver <- unres[,c(1,4)]

## Need to create polytomies. This function is too slow!!!
phy.poly <- to.make.poly(diver, tree.genus[[1]] )

to.make.poly <- function(diver, phy){
    ## Create a phylogeny with polytomies for clade data.
    ## diver is a matrix with first column is the tip.labels and second column is the total

    add.species <- function(tree, species, genus = NULL){
        ## This is just a simpler version of 'phytools' function add.species.to.genus.
        ## This function does not have the checkings. As a results is faster.
        ## Original 'phytools' function is more flexible.
    
        ## The 'where' parameter here is considered always == 'root'.
        ii <- grep(paste(genus, "_", sep = ""), tree$tip.label)
        if (length(ii) > 1) {
            nn <- findMRCA(tree, tree$tip.label[ii])
            tree <- bind.tip(tree, gsub(" ", "_", species), where = nn)
        } else {
            nn <- ii
            tree <- bind.tip(tree, gsub(" ", "_", species), where = nn, 
                             position = 0.5 * tree$edge.length[which(tree$edge[, 2] == nn)])
        }
        return(tree)
    }
    
    phy$tip.label <- paste(phy$tip.label,"_sp1",sep="")
    sp.list <- lapply(1:dim(diver)[1], FUN = function(x) paste(diver[x,1],"_sp",2:diver[x,2],sep=""))
    for(i in 1:length(sp.list)){
        print(paste("Clade ",i," of ",length(sp.list)," with ",diver[i,2]," species.",sep=""))
        for(j in 1:length(sp.list[[i]])){
            phy <- add.species(phy, species=sp.list[[i]][j], genus=diver[i,1])
        }
    }
    return(phy)
}
