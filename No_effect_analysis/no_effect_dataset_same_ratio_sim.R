## There is a major problem in the previous simulations.
## The datasets used a null model does not look like our data.
## The ratios are very off. The empirical proportion is 213 -> 1 and 381 -> 0.

## Here I am going to sample from a binomial distribuition with probability
##      proportional to the diversity of the trait across the tree.
## The result should be datasets in average more close to the real data and
##      also with not the same effect as the empirical data.

## Then I am going to use this distribuition of datasets across the empirical
##      trees and calculate the MLE estimate of the constrained and full
##      models.

library(diversitree)
source("../functions/analysis.R")

## Load BiSSE data:
load("../data/data_for_BiSSE.RData")
state <- as.numeric(st[,2])
names(state) <- st[,1]

## The observed probability of state 1:
prob <- 213 / (213 + 381)

## Use the 'to.make.rd.state' function to create randomly sampled datasets:
sim.st <- lapply(1:100, function(x) to.make.rd.state(state, unres, freq0 = 1 - prob, freq1 = prob) )

## Check the ratio of the simulated data:
## pp.unres.st1 <- lapply(1:100, function(x) sum(sim.st[[x]][[2]][,3]) )
## pp.st1 <- lapply(1:100, function(x) pp.unres.st1[[x]] + sum(sim.st[[x]][[1]], na.rm = TRUE) )
## pp.st1 <- do.call(cbind, pp.st1)[1,]

## hist(pp.st1 / (213 + 381), main = "States proportion", xlab = "St.1 / (St.0 + St.1)"
##    , breaks = 15)
## abline(v = prob, lty = 2)

save(sim.st, file = "Same_ratio_sim_data.RData")

## Start the analysis:
## Here I am running for just one tree. The tree number 1.

library(multicore)
index <- 1:100

sim.fit.cons <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tree.genus[[1]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , constrain = "TRUE", flag = paste("sim.fit.cons", x, sep=""))
                       , mc.cores = 20)
save(sim.fit.cons, file = "sim.fit.cons.RData")
rm(sim.fit.cons)

sim.fit.full <- mclapply(index, FUN = function(x) run.bisse.mle(tree = tree.genus[[1]]
                                  , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
                                  , constrain = "FALSE", flag = paste("sim.fit.full", x, sep=""))
                       , mc.cores = 20)
save(sim.fit.full, file = "sim.fit.full.RData")
rm(sim.fit.full)

#############################################
## Now get the MLE for the empirical dataset:
## Create state vector for BiSSE:

emp.fit.full <- run.bisse.mle(tree = tree.genus[[1]], st = state, unres = unres, constrain = "FALSE"
                            , flag = "Empirical - full")
save(emp.fit.full, file = "emp.fit.full.RData")

emp.fit.cons <- run.bisse.mle(tree = tree.genus[[1]], st = state, unres = unres, constrain = "TRUE"
                            , flag = "Empirical - cons")
save(emp.fit.cons, file = "emp.fit.cons.RData")
