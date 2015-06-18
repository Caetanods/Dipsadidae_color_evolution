## Script to calculate the Deviance Information Criteria and get results from model
##   selection.

require(diversitree)
require(parallel)

## Get the functions and the data:
source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")

## To calculate the DIC we need the likelihood for the parameters under both the contrained and
##     full models. Next steps will produce the likelihood function using the 'make.bisse'
##     function and calculate the likelihood given parameter estimates.

## Calculate for the trait-independent (simpler) BiSSE model:

foo.dic.one <- function(res, tree, st, unres){
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
    lik.c <- constrain(lik, lambda1~lambda0, mu1~mu0)
    b <- dim(res)[1]*50/100
    dic.one <- dic.mcmcsamples(res, burnin = b, lik = lik.c)
    return(dic.one)
}

## WARNING: This is a long process and the following line is set up to run using 15 cores.
dic.one <- mclapply(1:length(mcmc.onerate), FUN = function(x)
                  foo.dic.one(mcmc.onerate[[x]][[1]], tree.genus[[x]], st, unres), mc.cores = 15)

## Calculate for the trait-dependent (more complex) BiSSE model:

foo.dic.two <- function(res, tree, st, unres){
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
    b <- dim(res)[1]*50/100
    dic.two <- dic.mcmcsamples(res, burnin = b, lik = lik)
    return(dic.two)
}
dic.two <- mclapply(1:length(mcmc.onerate), FUN = function(x)
                  foo.dic.two(mcmc.onerate[[x]][[1]], tree.genus[[x]], st, unres), mc.cores = 15)

## Get the result as vectors:
dic.one <- as.vector( do.call(cbind, dic.one) )
dic.two <- as.vector( do.call(cbind, dic.two) )

save(dic.two, dic.one, file = "./data/dic_BiSSE_MCMC_results.RData")

##############################################################
## Analysis of the result of DIC and the model selection:

load("./data/dic_BiSSE_MCMC_results.RData")

dif <- dic.one - dic.two

## Quick graph to show the distribution of the DIC. (Not paired)
par(mfrow = c(1,2) )
hist(dic.one)
hist(dic.two)

range(dif)
mean(dif)

## The best model is the model with the smaller value of DIC.
## The following graph show the difference between the DIC of the simpler model - the complex model.
## If the complex model is the best, than values should be positive. The larger the difference,
##       the better the complex model is.

pdf("hist_DIC_results.pdf")
hist(dif, xlim = c(5,30), breaks = 20, freq = FALSE, main = "", xlab = "DIC difference values")
dev.off()
