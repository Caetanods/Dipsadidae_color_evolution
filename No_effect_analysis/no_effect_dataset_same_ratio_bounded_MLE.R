## This is a test to check if bounding the MLE estimate to parameters values that make sense
##      improve the MLE estimate and convergence.
## For this I will use the 'optim' optimizer with the bounded 'L-BFGS-B' method.

library(diversitree)
library(multicore)
source("../functions/analysis.R")

## Get sim data from previous run:
load("./Same_ratio_sim_data.RData")

## This is the list of runs which failed to converge in previous analysis.
## Will perform test to check if they get better.
fail <- c(1,3,4,6,8)

## Making the test with only the ones that failed.

sim.fit.cons <- mclapply(fail, FUN = function(x) run.bisse.mle(tree = tree.genus[[2]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , constrain = "TRUE", flag = paste("sim.fit.cons", x, sep=""))
                       , mc.cores = 20)
save(sim.fit.cons, file = "sim.fit.cons.bounded.RData")

sim.fit.full <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tree.genus[[2]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , constrain = "FALSE", flag = paste("sim.fit.full", x, sep=""))
                       , mc.cores = 20)
save(sim.fit.full, file = "sim.fit.full.bounded.RData")
