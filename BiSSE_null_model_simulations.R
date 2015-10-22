## This is similulations to test for the effect of the type I error on the model selection of the
##      BiSSE estimates. The idea is to simulated data where the effect of the traits on net diver-
##      sification is not present and test which is the best model. This analyses uses Likelihood
##      Ratio Tests in order to save (a lot!!) of computational resources.

## To generate the datasets we are sampling from a binomial distribution with probability of each
##      state proportional to the diversity of the trait in the empirical datasets.
## The result should be datasets in average similar to the real data and, since the data is
##      randomly distributed across the tips, with no effect of the traits in the diversification.

## Then I am going to use this distribuition of datasets across the empirical
##      trees and calculate the MLE estimate of the constrained and full
##      models. The BiSSE models fitted here are the same used to estimate the parameters under
##      the empirical analyses. The only difference is that a MLE method is used instead of a MCMC.
## Here we are simulating using 100 simulated datasets and across 10 trees, to account for phylogeny
##      uncertainty.

library(diversitree)
library(parallel)
source("./functions/analysis.R")

## Load BiSSE data:
load("./data/data_for_BiSSE.RData")
state <- as.numeric(st[,2])
names(state) <- st[,1]

## The observed probability of state 1:
prob <- 213 / (213 + 381)

## Use the 'to.make.rd.state' function to create randomly sampled datasets:
sim.st <- lapply(1:100, function(x) to.make.rd.state(state, unres, freq0 = 1 - prob, freq1 = prob) )

## Start the analysis:
## Declare function to run simulations for one tree.
## This function depend on the objects already loaded/created in the working space.
## This function will not work in another environment.
## WARNING: This simulation WILL OPEN NEW DIRECTORIES in the working directory. Please proceed
##      with caution.
## Uncomment the following block to run the analysis. Otherwise skip and load results.

## run.lrt.sim <- function(ii){
##     ## ii; index of the phylogeny to be used.
##     ## All other objects are already in the workspace.
##     ## function return void. all objects are saved to the directory.

##     dd <- paste("./tree_genus_", ii, sep = "")
##     dir.create( dd )
    
##     sim.fit.cons <- mclapply(1:100, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[ii]]
##                                   , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
##                                   , low = 0, up = 20
##                                   , constrain = "TRUE", flag = paste("sim.fit.cons.", x,"phy.",ii, sep=""))
##                        , mc.cores = 22)
##     save(sim.fit.cons, file = paste(dd, "/", "sim.fit.cons.RData", sep = "") )
##     rm(sim.fit.cons)
    
##     sim.fit.full <- mclapply(1:100, FUN = function(x) run.bisse.bound.mle(tree = tree.genus[[ii]]
##                                   , st = sim.st[[x]][[1]], unres = sim.st[[x]][[2]]
##                                   , low = 0, up = 20
##                                   , constrain = "FALSE", flag = paste("sim.fit.full", x,"phy.",ii, sep=""))
##                        , mc.cores = 22)
##     save(sim.fit.full, file = paste(dd, "/", "sim.fit.full.RData", sep = "") )
##     rm(sim.fit.full)

##     emp.fit.full <- run.bisse.bound.mle(tree = tree.genus[[ii]], st = state, unres = unres, constrain = "FALSE"
##                                       , low = 0, up = 20
##                                       , flag = paste("emp.fit.full.phy.", ii, sep="") )
##     save(emp.fit.full, file = paste(dd, "/", "emp.fit.full.RData", sep = "") )
##     rm(emp.fit.full)
    
##     emp.fit.cons <- run.bisse.bound.mle(tree = tree.genus[[ii]], st = state, unres = unres, constrain = "TRUE"
##                                       , low = 0, up = 20
##                                       , flag = paste("emp.fit.cons.phy.", ii, sep="") )
##     save(emp.fit.cons, file = paste(dd, "/", "emp.fit.cons.RData", sep = "") )
##     rm(emp.fit.cons)
## }

## sapply(2:10, function(x) run.lrt.sim(x) )

###########################################################################################
## SIMULATION RESULTS

## The following block will process and read the results of the simulation.
## You will need to uncomment and run the block if you completed the simulations above.
## Otherwise, skip the block and load the results.

## Get results from the different runs:
## dd <- paste("./no-effect-results/tree_genus_", 1:10, "/", sep="")
## res <- list()

## for(i in dd){
##     dat <- dir(i)
##     for(j in dat) load( paste(i, j, sep = "") )
##     ll <- list(emp.fit.cons, emp.fit.full, sim.fit.cons, sim.fit.full)
##     rm(emp.fit.cons, emp.fit.full, sim.fit.cons, sim.fit.full)
##     res[[i]] <- ll
## }

## ## This is the structure of the res object:
## ## For each j tree and k (1:100) simulation we have:
## #res[[j]][[1]] = emp.fit.cons
## #res[[j]][[2]] = emp.fit.full
## #res[[j]][[3]][[k]] = sim.fit.cons
## #res[[j]][[4]][[k]] = sim.fit.full

## ## Make the LRT:
## sim.lrt <- list()
## emp.lrt <- list()
## p.sim <- list()
## p.emp <- rep(NA, times = 10)

## for( i in 1:10 ){ ## Loop over trees.
##     sim.lrt[[i]] <- lapply(1:100, function(x) anova(res[[i]][[4]][[x]]$MLE, cons = res[[i]][[3]][[x]]$MLE ) )
##     emp.lrt[[i]] <- anova(res[[i]][[2]]$MLE, cons = res[[i]][[1]]$MLE )
##     p.sim[[i]] <- sapply(1:100, function(x) sim.lrt[[i]][[x]][[5]][2] )
##     p.emp[i] <- emp.lrt[[i]][[5]][2]
## }

## Save results:
## save(emp.lrt, sim.lrt, p.emp, p.sim, file = "no_effect_results.RData")

## Load the results from the simulations:
load("./data/no_effect_results.RData")

## TO CHECK THE RESULTS FROM THE ALTERNATIVE MODELS WE JUST NEED SOME DISTRIBUTION OF P VALUES.
## THE SIMULATION RESULTS CAN BE THE SAME.

## We are only interested in the case in which the p value for the empiric LRT is LOWER than
##    the 95% CI of the p values simulated under the null model. For this I am going to check
##    the 0.05 quantile of the distribution. This is equivalent to a single-tail Monte Carlo test.
qq.sim <- sapply(1:10, function(x) quantile(p.sim[[x]], probs = 0.05) )
check <- data.frame(p.sig = qq.sim, p.emp = p.emp )

## Look the results:
check
## Which empirical p value is LOWER than the simulation threshold?
check$p.sig > check$p.emp

## Make some plots for visualization of the data.
## Here the red line is the empirical p value from the LRT using the real dataset.
## The green line is p = 0.05, the most used p value threshould in ecology and evolution studies.
par(mfrow = c(3,2))
for( i in 1:5 ){
    hist(p.sim[[i]], breaks = 15, col = "grey", border = "white")
    abline(v = p.emp[i], col = "red", lwd = 2)
    abline(v = 0.05, col = "green", lty = 3, lwd = 2)
}
