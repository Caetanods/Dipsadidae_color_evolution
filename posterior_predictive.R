## This is the posterior predictive simulations for BiSSE analysis.

## Here I make draws of the joint posterior distribution of paramater values (each
##    draw is a given generation of the posterior).
## Then I create a vector of parameters from these draws.
## I simulate datasets under this parameters for the full model of BiSSE.
## I use summary statistics to compare the simulated datasets with the posterior distributions.
## Need to explore more summary statistics. (Maybe).

library(diversitree)

## Read results. Get combined posterior for full model with 50% of burnin.

## Uncomment following lines to download result file.
## Download the results file from "http://files.figshare.com/1696849/results_100_phylo_bisse.RData"
##    and use load() to get it to R.
#download.file("http://files.figshare.com/1696849/results_100_phylo_bisse.RData", "result_data.RData")
#load("result_data.RData")

load("./results/results_100_phylo_bisse.RData")

cc <- lapply(1:100, FUN = function(x) mcmc.tworate[[x]][[1]][5000:10000,c(-1,-8)])
comb.post <- as.matrix(do.call(rbind, cc), rownames.force = FALSE)

## Remove extra objects:
rm("mcmc.tworate","mcmc.onerate")

## Function to do the posterior predictive simulation:
post.predict.bisse <- function(post, max.taxa, size){
    ## Post is a matrix with the BiSSE paramters of the posterior.
    ## Size is the number of generations to be sampled from post.
    gen <- sample(1:nrow(post), size = size)    
    make.sim <- function(pars, max.taxa){
        repeat{
            phy <- tree.bisse(pars, max.taxa)
            if(length(phy$tip.label) == max.taxa){break}
        }
        ## Calculate series of summary stats.
        tl <- unname(diag(vcv.phylo(phy))[1])
        p <- sum(phy$tip.state) / length(phy$tip.state)
        vec <- c(p, tl)
        return(list(phy, vec))
    }
    res <- lapply(gen, FUN = function(x) make.sim(post[x,], max.taxa) )
    return(res)
}

## How many datasets to be simulated?
taxa <- 594
size <- 100

sim <- post.predict.bisse(post = comb.post, max.taxa = taxa, size = size)

## It works great. The questions is: What summary stats should I use?
sum.stat <- lapply(1:length(sim), function(x) sim[[x]][[2]] )
sum.stat <- do.call(rbind, sum.stat)
colnames(sum.stat) <- c("Freq_state_1","Tree_depth")
head(sum.stat)

## Import the data necessary to compare the summary stats.
load("./data/data_for_BiSSE.RData")
prop <- sum(unres$n1) / ( sum(unres$n0) + sum(unres$n1) )

## What other summary stats could be useful?
## 1) Frequency of state 1 (or 0) across the clades.
## 2) Number of changes from 0 to 1 and 0 to 1 using stochastic map.
## 3) Time (branch length) spent in 1 (or 0) from stochastic map.
## 4) Mean waiting times for speciation while in state 0 and state 1. For this
##    we need to use the stochastic map. We need to exclude the tips, since
##    they are still waiting to speciate. Then we paint the tree for branches
##    with the 0 and 1 states. Then we sum the waiting times. Waiting time is
##    the branch length from transition to state X until the next speciation
##    event. Compare the mean waiting time for state 1 and state 0.
