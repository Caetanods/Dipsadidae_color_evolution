## Here I am checking the results of the simulations where the tips were draw from a
##    distribution with equal ratio between the traits to the observed states.
## Then I calculate the log lik ratio test for them and compare to the empirical
##    lik ratio.

library(diversitree)
source("../functions/analysis.R")

## Get results:
ll.dir <- dir("./no-effect-results/MLE_shuffle_equal_ratio_emp", pattern = "*.RData")
ll.dat <- paste("./no-effect-results/MLE_shuffle_equal_ratio_emp/", ll.dir, sep = "")
for(i in ll.dat) load(i)

## Compare the models in the simulated data:
## Some of the analysis did not converged.
## The MLE that converged are c(2,5,7,9,10)
cc <- c(2,5,7,9,10)
av <- lapply(cc, function(x) anova( sim.fit.full[[x]]$MLE, cons=sim.fit.cons[[x]]$MLE ) )
## Get delta AIC; cons AIC - full AIC; larger values prefer full model.
delta.aic <- sapply(1:length(av), function(x) av[[x]][3][2,] - av[[x]][3][1,] )
delta.aic

## Make a lik ratio test of the empirical data:
av.emp <- anova(emp.fit.full$MLE, cons=emp.fit.cons$MLE)

## delta AIC:
av.emp[3][2,] - av.emp[3][1,]

## Making the simulation only 5x I can get a difference as strong or more strong than the
##      empirical. If I extend this simulations I am sure my empirical result
##      will fall within the expected for a random dataset.

