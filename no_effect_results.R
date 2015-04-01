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

###########################################
## Now we are going to compare the results with MCMC posterior distributions for
##        the 20 first set of trees.
## load the results:

load("no-effect-results/no_effect_results.RData")
load("no-effect-results/no_effect_sim_states.RData")

mcmc_post <- lapply(mcmc_res, function(x) x[5000:10000,])
net0 <- lapply(mcmc_post, function(x) x[,3] - x[,5] )
net1 <- lapply(mcmc_post, function(x) x[,4] - x[,6] )

## Make plot of the net diversification results from the simulations.
## Red is the Contrasting state and blue is the cryptic.
pdf("no_effect_results.pdf", width = 14, height = 14)
par(mfrow = c(4,5))
for(i in 1:length(net0) ){
    mn <- min( c(net0[[i]],net1[[i]]) )
    mx <- max( c(net0[[i]],net1[[i]]) )
    plot(density(net0[[i]], from = mn, to = mx), xlim = c(mn,mx), col = "blue"
       , main = paste("No effect sim #",i,sep = ""), xlab = "Net diversification", lwd = 1.5)
    lines(density(net1[[i]], from = mn, to = mx), col = "red", lwd = 1.5)
    abline(v = mean(net0[[i]]), col = "blue", lwd = 1.5, lty = 3)
    abline(v = mean(net1[[i]]), col = "red", lwd = 1.5, lty = 3)
}
dev.off()

## How big is the effect?
## Plot the distribution of the difference of net diversification in the simulated data
##      and in the empirical results.

## OBS: Getting only first 20 trees from the empirical data.
emp.20 <- mcmc.tworate[1:20]
emp.20 <- lapply(emp.20, function(x) x[[1]][5000:10000,] )

sim.dif <- lapply(mcmc_post, function(x) (x[,4] - x[,6]) - (x[,3] - x[,5]) )
sim.dif <- c(do.call(rbind, sim.dif))
emp.dif <- lapply(emp.20, function(x) (x[,3] - x[,5]) - (x[,2] - x[,4]) )
emp.dif <- c(do.call(rbind, emp.dif))

mn <- min(c(sim.dif, sim.dif))
mx <- max(c(sim.dif, sim.dif))

pdf("effect_difference.pdf")
plot(density(emp.dif, from = mn, to = mx), xlim = c(mn,mx), col = "red"
     , xlab = "Net diversification difference (st.1 - st.0)", lwd = 1.5
     , main = "")
hist(emp.dif, xlim = c(mn,mx), border = "red", xlab = "", main = "", add = TRUE, freq = FALSE)
lines(density(sim.dif, from = mn, to = mx), col = "blue", lwd = 1.5)
hist(sim.dif, xlim = c(mn,mx), border = "blue", xlab = "", main = "", add = TRUE, freq = FALSE)
abline(v = mean(emp.dif), col = "red", lwd = 1.5, lty = 3)
abline(v = mean(sim.dif), col = "blue", lwd = 1.5, lty = 3)
legend(-10, 0.4, legend = c("Simulation","Empirical"), lwd = 1.5, col = c("blue","red") )
dev.off()

## Are the distributions different?
## Alternative, x is shifted to the right of y.
wilcox.test(emp.dif, sim.dif, alternative = "greater")
## Not sure if I can test things using a test like that.

## Check the predictive concordance ( Gelfand, 1996 ).
outlier <- quantile(sim.dif, probs = 1-0.025)
sum(emp.dif <= outlier) / length(emp.dif)
get.overlap(emp.dif, sim.dif)

## The 0.98 proportion of the empirical distribution overlaps with the
##     simulated distribution. This means that the posterior distributions
##     of differences overlaps. The magnitude of the effect found in the empirical
##     data is within the expected if there is no effect in the data.

## How about the values of the parameters?
## How much overlap there is between the net diversification of the state 1 and state 0?
sim.net0 <- c(do.call(rbind, net0))
sim.net1 <- c(do.call(rbind, net1))
emp.net0 <- lapply(emp.20, function(x) x[,2] - x[,4] )
emp.net0 <- c(do.call(rbind, emp.net0))
emp.net1 <- lapply(emp.20, function(x) x[,3] - x[,5] )
emp.net1 <- c(do.call(rbind, emp.net1))

get.overlap(emp.net0, sim.net0)
get.overlap(emp.net1, sim.net1)

quantile(sim.net0, 0.975)
quantile(sim.net0, 0.025)

## This function does not work yet!!
get.overlap <- function(x,y){
    ## Calculate the portion of x that is inside the 95% HDI of y.
    ## This is the Predictive Concordance ( Gelfand, 1996 ).
    right <- x[which(x >= quantile(y, probs = 0.025))]
    left <- right[which(right <= quantile(y, probs = 0.095))]
    res <- length(left) / length(x)
    return(res)
}
