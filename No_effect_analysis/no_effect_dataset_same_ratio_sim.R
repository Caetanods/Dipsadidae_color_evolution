## Here I am going to sample from a binomial distribuition with probability
##      proportional to the diversity of the trait across the tree.
## The result should be datasets in average more close to the real data and
##      also with not the same effect as the empirical data.

## Then I am going to use this distribuition of datasets across the empirical
##      trees and calculate the MLE estimate of the constrained and full
##      models.

## I am making parameter estimates for 100 simulated datasets for each empirical tree and 10 empirical trees.

library(diversitree)
library(multicore)
source("../functions/analysis.R")

## Load BiSSE data:
load("../data/data_for_BiSSE.RData")
state <- as.numeric(st[,2])
names(state) <- st[,1]

## The observed probability of state 1:
prob <- 213 / (213 + 381)

## Use the 'to.make.rd.state' function to create randomly sampled datasets:
sim.st <- lapply(1:100, function(x) to.make.rd.state(state, unres, freq0 = 1 - prob, freq1 = prob) )
save(sim.st, file = "Same_ratio_sim_data.RData")

## Start the analysis:

## Declare function to run simulations for one tree:

run.lrt.sim <- function(ii){
    ## ii; index of the phylogeny to be used.
    ## All other objects are already in the workspace.
    ## function return void. all objects are saved to the directory.

    dd <- paste("./tree_genus_", ii, sep = "")
    dir.create( dd )
    
    sim.fit.cons <- mclapply(1:100, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[ii]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , low = 0, up = 20
                                  , constrain = "TRUE", flag = paste("sim.fit.cons.", x,"phy.",ii, sep=""))
                       , mc.cores = 22)
    save(sim.fit.cons, file = paste(dd, "/", "sim.fit.cons.RData", sep = "") )
    rm(sim.fit.cons)
    
    sim.fit.full <- mclapply(1:100, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[ii]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , low = 0, up = 20
                                  , constrain = "FALSE", flag = paste("sim.fit.full", x,"phy.",ii, sep=""))
                       , mc.cores = 22)
    save(sim.fit.full, file = paste(dd, "/", "sim.fit.full.RData", sep = "") )
    rm(sim.fit.full)

    emp.fit.full <- run.bisse.bound.mle(tree = tree.genus[[ii]], st = state, unres = unres, constrain = "FALSE"
                                      , low = 0, up = 20
                                      , flag = paste("emp.fit.full.phy.", ii, sep="") )
    save(emp.fit.full, file = paste(dd, "/", "emp.fit.full.RData", sep = "") )
    rm(emp.fit.full)
    
    emp.fit.cons <- run.bisse.bound.mle(tree = tree.genus[[ii]], st = state, unres = unres, constrain = "TRUE"
                                      , low = 0, up = 20
                                      , flag = paste("emp.fit.cons.phy.", ii, sep="") )
    save(emp.fit.cons, file = paste(dd, "/", "emp.fit.cons.RData", sep = "") )
    rm(emp.fit.cons)
}

sapply(2:10, function(x) run.lrt.sim(x) )

## The tree number 1 is already done in macumba in the following directory:
## '/Documents/snake_manuscript/diversitree_analysis/mle_equal_prop_tree_gen_1'
