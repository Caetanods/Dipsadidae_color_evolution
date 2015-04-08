## This is a test to check if bounding the MLE estimate to parameters values that make sense
##      improve the MLE estimate and convergence.
## For this I will use the 'optim' optimizer with the bounded 'L-BFGS-B' method.

library(diversitree)
library(multicore)
source("../functions/analysis.R")

## Load BiSSE data:
load("../data/data_for_BiSSE.RData")
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Get sim data from previous run:
load("./Same_ratio_sim_data.RData")

## This is the list of runs which failed to converge in previous analysis.
## Will perform test to check if they get better.
fail <- c(1,3,4,6,8)

## Just a test of the bounded MLE:
## test <- run.bisse.bound.mle(tree = tree.genus[[1]], st = sim.st[[1]][[1]], low = 0, up = 10
##                   , unres = sim.st[[1]][[2]], constrain = "TRUE"
##                   , flag = paste("sim.fit.cons", 1, sep=""))

## Making the test with only the ones that failed.

tasks <- list(
    sim.fit.cons <- function() mclapply(fail, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[1]]
                                  , st = sim.st[[x]][[1]], low = 0, up = 10
				  , unres = sim.st[[x]][[2]], constrain = "TRUE"
                                  , flag = paste("sim.fit.cons", x, sep=""))
                       , mc.cores = 5),
    sim.fit.full <- function() mclapply(fail, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[1]]
                                  , st = sim.st[[x]][[1]], low = 0, up = 10
                                  , unres = sim.st[[x]][[2]], constrain = "FALSE"
                                  , flag = paste("sim.fit.full", x, sep=""))
                       , mc.cores = 5)
)

out <- mclapply(tasks, function(f) f(), mc.cores = 20)

save(out, file = "same.ratio.bound.mle.test.RData")
