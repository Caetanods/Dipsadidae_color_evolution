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

## Make the LRT:
lrt.sim.free <- lapply(cc, function(x) anova( sim.fit.full[[x]]$MLE, cons=sim.fit.cons[[x]]$MLE ) )
p.sim.free <- sapply(1:length(lrt.sim.free), function(x) lrt.sim.free[[x]][[5]][2] )
lrt.emp <- anova(emp.fit.full$MLE, cons=emp.fit.cons$MLE)
p.emp <- lrt.emp[[5]][2]

## Making the simulation only 5x I can get a difference as strong or more strong than the
##      empirical. If I extend this simulations I am sure my empirical result
##      will fall within the expected for a random dataset.

## I bounded the MLE estimate and re-estimated the parameters using the simulations
##      which did not converged.
sim.fit.cons.bound <- out[[1]]
sim.fit.full.bound <- out[[2]]
ll <- length(sim.fit.cons.bound)
lrt.sim.bound <- lapply(1:ll, function(x) anova(sim.fit.full.bound[[x]]$MLE, cons=sim.fit.cons.bound[[x]]$MLE))
p.sim.bound <- sapply(1:ll, function(x) lrt.sim.bound[[x]][[5]][2] )

c(p.sim.bound, p.sim.free) > p.emp
hist( c(p.sim.bound, p.sim.free) )
abline(v = p.emp)

pars <- sapply(1:ll, function(x) coef(sim.fit.full.bound[[x]]$MLE) )
pars
