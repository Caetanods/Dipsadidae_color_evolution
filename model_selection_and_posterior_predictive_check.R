## Script to calculate the Deviance Information Criteria for model selection and Posterior
##   predictive checks for such models.

################################################################################################
################################################################################################
## DEVIANCE INFORMATION CRITERIA:
##   Model selection here is performed for all the different categorizations.
##   Note that results from the DIC model selection do not include the correction for the
##           bias of the BiSSE model in selecting the more complex model even when the correct
##           model is the simpler one. Simulations to correct for this bias were made in other
##           scripts.

require(diversitree)
require(parallel)

## Get the functions and the data:
source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")
## The dataset bellow need to be downloaded from FigShare first. If you do not have the RData. Then:
## download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData"
##             , destfile = "./data/results_100_phylo_bisse.RData", method = "wget")
load("./data/results_100_phylo_bisse.RData")

## To calculate the DIC we need the likelihood for the parameters under both the contrained and
##     full models. Next steps will produce the likelihood function using the 'make.bisse'
##     function and calculate the likelihood given parameter estimates.
## Since the DIC involves estimating the likelihood for the parameters again. The estimation takes
##     some time. Uncomment the block bellow to run the complete analysis again. Otherwise, skip
##     the block to load the result object.

## WARNING: This is a long process and the following line is set up to run using 15 cores.
## Uncomment to run.

## ## Calculate DIC for the trait-independent (simpler) BiSSE model:
## foo.dic.one <- function(res, tree, st, unres){
##     tree <- multi2di(tree)
##     mmm <- match(tree$tip.label, names(st))
##     lik <- make.bisse(tree, st[mmm], unresolved = unres)
##     lik.c <- constrain(lik, lambda1~lambda0, mu1~mu0)
##     b <- dim(res)[1]*50/100
##     dic.one <- dic.mcmcsamples(res, burnin = b, lik = lik.c)
##     return(dic.one)
## }

## dic.one <- mclapply(1:length(mcmc.onerate), FUN = function(x)
##                   foo.dic.one(mcmc.onerate[[x]][[1]], tree.genus[[x]], st, unres), mc.cores = 15)

## ## Calculate for the trait-dependent (more complex) BiSSE model:

## foo.dic.two <- function(res, tree, st, unres){
##     tree <- multi2di(tree)
##     mmm <- match(tree$tip.label, names(st))
##     lik <- make.bisse(tree, st[mmm], unresolved = unres)
##     b <- dim(res)[1]*50/100
##     dic.two <- dic.mcmcsamples(res, burnin = b, lik = lik)
##     return(dic.two)
## }
## dic.two <- mclapply(1:length(mcmc.onerate), FUN = function(x)
##                   foo.dic.two(mcmc.onerate[[x]][[1]], tree.genus[[x]], st, unres), mc.cores = 15)

## ## Get the result as vectors:
## dic.one <- as.vector( do.call(cbind, dic.one) )
## dic.two <- as.vector( do.call(cbind, dic.two) )

## save(dic.two, dic.one, file = "./data/dic_BiSSE_MCMC_results.RData")

##############################################################
## Analysis of the result of DIC and the model selection:

## Load the results of the DIC estimate.
load("./data/dic_BiSSE_MCMC_results.RData")

## Since the results for the DIC are produced from analysis over 100 phylogenies from the posterior
##     we can make a paired comparison of the DIC scores.
dif <- dic.one - dic.two

## Quick graph to show the distribution of the differences in the DIC scores.
hist(dif, probability = TRUE, main="", xlab = "DIC difference", col = "grey", border = "white")
lines(density(dif), col = "red", lwd = 1.5)
abline(v=mean(dif), col = "blue", lty = 4, lwd = 2)
abline(v=median(dif), col = "green", lty = 4, lwd = 2)

## The best model is the model with the smaller value of DIC.
## If the complex model is the best, than values should be positive. The larger the difference,
##       the better the complex model is.
## The above graph show that the full BiSSE model, where the net diversification is dependent of the
##       color pattern of the lineages is the best model.
## Further analysis will show that this difference is a result of the bias of the BiSSE model.
## Even if data is produced under a model with no difference between the net diversification associa
##       ted with each of the color states, we will recover the best model as the more parameter
##       rich, the full BiSSE model.

############################################################################
## DIC ANALYSIS AND RESULTS FOR THE ALTERNATIVE CATEGORIES:


################################################################################################
################################################################################################
## POSTERIOR PREDICTIVE CHECKS

## The following results assumes the main categorization, which bright and contrasting colors are
##     set as the contrasting (1) states, including species in which the colors are in the venter
##     and in case that ontogenetic changes occur, the data come from the juveniles. All species
##     that do not show contrasting colors are defined as the cryptic state (0).

## If the following RData is not in the "./data" directory. Download the file from FigShare:
## download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
load("./data/results_100_phylo_bisse.RData")

## Get parameters estimates from the full model and constrained models. Take out 50% burnin.
two.rate.mean <- list()
for(i in 1:100){
    two.rate.mean[[i]] <- mcmc.tworate[[i]][[1]][5000:10000,]
}
two.rate.mean <- do.call(rbind, two.rate.mean)[,c(-1,-8)]

one.rate.mean <- list()
for(i in 1:100){
    one.rate.mean[[i]] <- mcmc.onerate[[i]][[1]][5000:10000,]
}
one.rate.mean <- do.call(rbind, one.rate.mean)[,c(-1,-6)]

## Load results of the MLE ancestral state estimate under the BiSSE model:
## For more details see 'make_figures.R'.
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

## Source functions:
source("./functions/analysis.R")

## WARNING: The analysis bellow can take some time. Uncomment the block to run the analysis.
## Otherwise, skip the block and load the results.

## ## Simulating data under the full model. Parameters are sampled from the joint posterior distribution
## ##    across 100 sampled trees. Tree depth set to 1.0, equal to empirical tree.
## ## The function 'sims' makes the simulations. Check './functions/analysis.R' for more info.

## sim.full <- lapply(1:1000, function(x) sims(pars = two.rate.mean, time = 1.0, freq.root = mean.root.two.rate) )
## sim.full.st0 <- sapply(1:1000, function(x)
##     (length(sim.full[[x]]$tip.state) - sum(sim.full[[x]]$tip.state)) / length(sim.full[[x]]$tip.state)
##                        )
## sim.full.rich <- sapply(1:1000, function(x) length(sim.full[[x]]$tip.state) )

## ## Making the same set of simulations for the constrained model:
## ## Need to create a vector of paramters with six elements.
## ## The order of the parameters is lambda0; lambda1; mu0; mu1; q01; q10.
## one.rate.par <- cbind( one.rate.mean[,1], one.rate.mean[,1]
##                      , one.rate.mean[,2], one.rate.mean[,2]
##                      , one.rate.mean[,3:4]
##                       )
## sim.con <- lapply(1:1000, function(x) sims(pars = one.rate.par, time = 1.0, freq.root = mean.root.one.rate) )
## sim.con.st0 <- sapply(1:1000, function(x)
##     (length(sim.con[[x]]$tip.state) - sum(sim.con[[x]]$tip.state)) / length(sim.con[[x]]$tip.state)
##                        )
## sim.con.rich <- sapply(1:1000, function(x) length(sim.con[[x]]$tip.state) )

## ## Save results.
## save(sim.full.st0, sim.full.rich, sim.con.st0, sim.con.rich, file = "./data/post_check.RData")

## Load the results:

load("./data/post_check.RData")

##########################################################################
## Plotting the results:

par(mfrow = c(2,2))
hist(sim.full.rich, breaks = 25, main = "Total richness - Full model", xlab = "Number of species"
    , xlim = c(0,2000), col = "grey", border = "white")
abline(v = 594, col = "red", lty = 2, lwd = 1.5)
legend(x = 1250, y = 450, legend = "Observed", col = "red", lty = 2, lwd = 1.5, bty = "n")
hist(sim.full.st0, breaks = 25, main = "State frequency - Full model", xlab = "Frequency of state 0"
    , xlim = c(0,1), col = "grey", border = "white" )
abline(v = 0.64, col = "blue", lty = 2, lwd = 1.5)
legend(x = 0, y = 140, legend = "Observed", col = "blue", lty = 2, lwd = 1.5, bty = "n")

hist(sim.con.rich, breaks = 25, main = "Total richness - Constrained model", xlab = "Number of species"
   , xlim = c(0,2000), col = "grey", border = "white" )
abline(v = 594, col = "red", lty = 2, lwd = 1.5)
legend(x = 1250, y = 200, legend = "Observed", col = "red", lty = 2, lwd = 1.5, bty = "n")
hist(sim.con.st0, breaks = 25, main = "State frequency - Constrained model", xlab = "Frequency of state 0"
    , xlim = c(0,1), col = "grey", border = "white" )
abline(v = 0.64, col = "blue", lty = 2, lwd = 1.5)
legend(x = 0, y = 140, legend = "Observed", col = "blue", lty = 2, lwd = 1.5, bty = "n")

## What is the proportion of the distribution more extreme than the observed values?
## Full model:
length( which(sim.full.rich > 594) ) / 1000
length( which(sim.full.st0 > 0.64) ) / 1000
## Constrained model:
length( which(sim.con.rich > 594) ) / 1000
length( which(sim.con.st0 > 0.64) ) / 1000

################################################################################################
################################################################################################
## POSTERIOR PREDICTIVE CHECKS FOR THE ALTERNATIVE CATEGORIZATIONS

## The posterior predictive analyses for the other 3 alternative categorizations.
