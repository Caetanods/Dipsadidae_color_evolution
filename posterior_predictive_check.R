## Perform the posterior predictive checks for the BiSSE analysis under the constrained and
##    unconstrained MCMC parameter estimates.

library(diversitree)

## Load BiSSE parameter estimates.
## Download the file from FigShare.
## download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
## load("./results_bisse.RData")
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")

## Get parameters estimates from the full model and constrained models. Take out 50% burnin.

two.rate.mean <- list()
for(i in 1:100){
    two.rate.mean[[i]] <- apply(mcmc.tworate[[i]][[1]][5000:10000,], 2, mean)
}
two.rate.mean <- do.call(rbind, two.rate.mean)[,c(-1,-8)]
two.rate.mean <- apply(two.rate.mean, 2, mean)

one.rate.mean <- list()
for(i in 1:100){
    one.rate.mean[[i]] <- apply(mcmc.onerate[[i]][[1]][5000:10000,], 2, mean)
}
one.rate.mean <- do.call(rbind, one.rate.mean)[,c(-1,-6)]
one.rate.mean <- apply(one.rate.mean, 2, mean)

## Load results of the MLE ancestral state estimate under the BiSSE model:
load("./data/asr_100_phylo.RData")

## Get the distribution of ASR of the root for all the trees:
root.two.rate <- t(sapply(asr.two.rate, FUN = function(x) x[,1]))
colnames(root.two.rate) <- c("State0","State1")
mean.root.two.rate <- apply(root.two.rate, 2, mean)

root.one.rate <- t(sapply(asr.one.rate, FUN = function(x) x[,1]))
colnames(root.one.rate) <- c("State0","State1")
mean.root.one.rate <- apply(root.one.rate, 2, mean)

##############################################################
## Posterior predictive simulations.
## Define functions:

sims <- function(pars, time, freq.root){
    ## Function that always return a tree with given parameters.
    ## This function will stop the BiSSE tree simulation using time only.
    ## The root state is drawn for a binomial distribution.
    repeat
        {
            root <- sample(c(0,1), size = 1, prob = freq.root)
            phy <- tree.bisse(pars, max.t = time, x0 = root)
            if(!is.null(phy))
                {
                    break
                }
        }
    return(phy)
}

## Simulating data under the full model. Parameters equal to mean parameter estimates
##   across the 100 sampled trees. Tree depth set to 1.0, equal to empirical tree.
## The order of the parameters is lambda0; lambda1; mu0; mu1; q01; q10.

sim.full <- lapply(1:1000, function(x) sims(pars = two.rate.mean, time = 1.0, freq.root = mean.root.two.rate) )
sim.full.st0 <- sapply(1:1000, function(x)
    (length(sim.full[[x]]$tip.state) - sum(sim.full[[x]]$tip.state)) / length(sim.full[[x]]$tip.state)
                       )
sim.full.rich <- sapply(1:1000, function(x) length(sim.full[[x]]$tip.state) )

## Making the same set of simulations for the constrained model:
## Need to create a vector of paramters with six elements.
## The order of the parameters is lambda0; lambda1; mu0; mu1; q01; q10.
one.rate.par <- c(one.rate.mean[1],one.rate.mean[1],one.rate.mean[2], one.rate.mean[2], one.rate.mean[3:4])
sim.con <- lapply(1:1000, function(x) sims(pars = one.rate.par, time = 1.0, freq.root = mean.root.one.rate) )
sim.con.st0 <- sapply(1:1000, function(x)
    (length(sim.con[[x]]$tip.state) - sum(sim.con[[x]]$tip.state)) / length(sim.con[[x]]$tip.state)
                       )
sim.con.rich <- sapply(1:1000, function(x) length(sim.con[[x]]$tip.state) )

##########################################################################
## Plotting the results:

pdf("Posterior_check.pdf", width = 9, height = 9)

par(mfrow = c(2,2))
hist(sim.full.rich, breaks = 25, main = "Total richness - Full model", xlab = "Number of species"
    , xlim = c(0,2000) )
abline(v = 594, col = "red", lty = 2, lwd = 1.5)
legend(x = 1250, y = 450, legend = "Observed", col = "red", lty = 2, lwd = 1.5, bty = "n")
hist(sim.full.st0, breaks = 25, main = "State frequency - Full model", xlab = "Frequency of state 0"
    , xlim = c(0,1) )
abline(v = 0.64, col = "blue", lty = 2, lwd = 1.5)
legend(x = 0, y = 140, legend = "Observed", col = "blue", lty = 2, lwd = 1.5, bty = "n")

hist(sim.con.rich, breaks = 25, main = "Total richness - Constrained model", xlab = "Number of species"
   , xlim = c(0,2000) )
abline(v = 594, col = "red", lty = 2, lwd = 1.5)
legend(x = 1250, y = 200, legend = "Observed", col = "red", lty = 2, lwd = 1.5, bty = "n")
hist(sim.con.st0, breaks = 25, main = "State frequency - Constrained model", xlab = "Frequency of state 0"
    , xlim = c(0,1) )
abline(v = 0.64, col = "blue", lty = 2, lwd = 1.5)
legend(x = 0, y = 140, legend = "Observed", col = "blue", lty = 2, lwd = 1.5, bty = "n")

dev.off()
