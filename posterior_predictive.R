## This is the posterior predictive simulations for BiSSE analysis.

## Here I make draws of the joint posterior distribution of paramater values (each
##    draw is a given generation of the posterior).
## Then I create a vector of parameters from these draws.
## I simulate datasets under this parameters for the full model of BiSSE.
## I use summary statistics to compare the simulated datasets with the posterior distributions.
## Need to explore more summary statistics. (Maybe).

library(diversitree)

## Read results. Get combined posterior for full model with 50% of burnin.
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
        p <- sum(phy$tip.state)
        vec <- c(p, tl)
        return(list(phy, vec))
    }
    res <- lapply(gen, FUN = function(x) make.sim(post[x,], max.taxa))
    return(res)
}

## How many datasets to be simulated?
taxa <- 594
size <- 10

sim <- post.predict.bisse(post = comb.post, max.taxa = taxa, size = size)

## It works great. The questions is: What summary stats should I use?
sim[[1]]

## Import the data necessary to compare the summary stats.
