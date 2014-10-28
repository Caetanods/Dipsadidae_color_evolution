## This script will do some simulations to check how the type 1 error affect my results.

## Sim1: Randomizing the data on the tips.
## For this I will create a new dataset of tip data generated with the same frequency of states of the
##     real dataset.

library(multicore)
library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Get total of species and frequency of each state.
## Sum of the 'unres' col plus the monotypic genera.
total.sp <- sum(unres$Nc) + length(which(!is.na(state)))
st1.freq <- ( sum(unres$n1) + length(which(state == "1")) ) / total.sp
st0.freq <- ( sum(unres$n0) + length(which(state == "0")) ) / total.sp

## Function to create the 'unres' and 'state' objects from sampled data:
rd.unres <- unres
rd.unres$n0 <- NA
rd.unres$n1 <- NA
init <- 1

for(i in 1:length(rd.unres$Nc)){
    nc <- rd.unres$Nc[i]
    rd <- sample(x = c(0,1), size = nc, replace = TRUE, prob = c(st0.freq, st1.freq))
    rd.state <- c(table(rd))
    print(rd.state)
    rd.unres$n0[i] <- as.numeric(rd.state[1])
    rd.unres$n1[i] <- as.numeric(rd.state[2])
    init <- init+nc
}

rd[1:sum(unres$Nc)]
head(unres)

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run example of analysis:
res.constrain <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Sim2: Simulating a no-effect MK2 model.
