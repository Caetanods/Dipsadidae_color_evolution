## Make analysis of the results of the no_effect_analysis.R

library(diversitree)
library(geiger)
source("./functions/data-prepare.R")
source("./functions/analysis.R")

## Load the results of simulations:
load("no-effect-results/no_effect_mle_results.RData")

###################################
## Check whether proportion of simulated states is representative of the empirical data.
## The empirical proportion is 213 -> 1 and 381 -> 0.

## Get proportion of states from sim datasets:
head(unres[[1]])
unres.st1 <- sapply(unres, function(x) sum(x[,3]) )
vector.st1 <- sapply(st.bisse, function(x) sum(x, na.rm = TRUE) )
st1 <- unres.st1 + vector.st1
prop.st1 <- st1 / (213 + 381)
prop.emp <- 213 / (213 + 381)

hist(prop.st1, main = "States proportion", xlab = "St.1 / (St.0 + St.1)", breaks = 15, xlim = c(0,1))
abline(v = prop.emp, lty = 2)
## The proportion of contrasting species in the simulated datasets is equal or lower than the observed value.
## This makes sense since the MK2 model is inadequade when there is a true effect of the state dependent
##      diversification.
## There is less contrasting species in the simulated datasets than what is observed, this means that
##      the higher speciation rates would compensante this difference, since when in state 1 lineages
##      accumulate faster.

###############################################
## Get the MLE of the parameter estimates to compare with the empirical results.

mle.no.effect[[1]][[2]]$par
mle.res <- lapply(mle.no.effect, function(x) x[[2]]$par)
mle.res <- do.call(rbind, mle.res)

net0 <- mle.res[,1] - mle.res[,3]
net1 <- mle.res[,2] - mle.res[,4]

## Get results from the empirical MCMC.
## To compare with the MLE we will use the mean of the posterior distributions
##     across the 100 trees.
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")
rm(mcmc.onerate)

mcmc.net0 <- sapply(1:100, function(x) mean( mcmc.tworate[[x]][[1]][,2] - mcmc.tworate[[x]][[1]][,4] ) )
mcmc.net1 <- sapply(1:100, function(x) mean( mcmc.tworate[[x]][[1]][,3] - mcmc.tworate[[x]][[1]][,5] ) )

## Plot the net diversification ML estimates across from the empirical and simulated results:
mn <- min( c(net0, net1, mcmc.net0, mcmc.net1) )
mx <- max( c(net0, net1, mcmc.net0, mcmc.net1) )
plot(density(net0, from = mn, to = mx), xlim = c(mn,mx), ylim = c(0,2.5), col = "grey"
   , xlab = "Net diversification (MLE)", lwd = 1.5, lty = 3
   , main = "Distribution of ML estimates")
lines(density(net1, from = mn, to = mx), col = "red", lwd = 1.5, lty = 3)
lines(density(mcmc.net0, from = mn, to = mx), col = "grey", lwd = 1.5)
lines(density(mcmc.net1, from = mn, to = mx), col = "red", lwd = 1.5)

## Plot the magnitude of the difference between the net diversifications from the empirical data
##      and the simulated.
mn <- min( c(net1 - net0, mcmc.net1 - mcmc.net0) )
mx <- max( c(net1 - net0, mcmc.net1 - mcmc.net0) )
plot(density(net1 - net0, from = mn, to = mx), xlim = c(mn,mx), ylim = c(0,1.5), col = "blue"
   , xlab = "Net diversification (MLE)", lwd = 1.5, main = "Magnitude of difference (net1 - net0)")
lines(density(mcmc.net1 - mcmc.net0, from = mn, to = mx), col = "red", lwd = 1.5)

#############################################################
#############################################################
## Check results of simulations using the transition rates from the contrained model
## This rates are more simmetrical than the full model estimates. Here we assume, when
##      estimating the transition rates between the states that there is no effect of
##      the different states in the rates of diversification.

## Clean workspace:
rm(list = ls())

## Get results:
load("./no-effect-results/MLE_constrained_model/no_effect_sim_data.RData")
load("./no-effect-results/MLE_constrained_model/no_effect_cons_mle_results.RData")

## Check proportion of simulated states. Empirical is 213 -> 1 and 381 -> 0.

cons.unres.st1 <- sapply(cons.unres, function(x) sum(x[,3]) )
cons.st1 <- cons.unres.st1 + sapply(cons.st.bisse, function(x) sum(x, na.rm = TRUE) )
cons.prop.st1 <- cons.st1 / (213 + 381)
prop.emp <- 213 / (213 + 381)

hist(cons.prop.st1, main = "States proportion", xlab = "St.1 / (St.0 + St.1)", breaks = 15, xlim = c(0,1))
abline(v = prop.emp, lty = 2)

########################################################
## The simulated datasets under both the constrained and the full models are very bad.
## The proportion of state 1 and 0 is very different of the proportion observed.
## We can try to get a more realistic simulation by estimating the MK model under the original
##    phylogeny in which all tips have branch lengths calculated under a uncorrelated BD model.
## For this I need to get the original trees from BEAST. Than estimate a simple ML2 model.
## Use the new estimate to simulate datasets and check if the propotion is more similar to the
##    empirical data.
## This is another alternative, not clear if this is in any means better than the past things I
##    tried.
