require(diversitree)
require(parallel)

source("../functions/analysis.R")

foo.dic.one <- function(res, tree, st, unres){
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
    lik.c <- constrain(lik, lambda1~lambda0, mu1~mu0)
    b <- dim(res)[1]*50/100
    dic.one <- dic.mcmcsamples(res, burnin = b, lik = lik.c)
    return(dic.one)
}

foo.dic.two <- function(res, tree, st, unres){
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
    b <- dim(res)[1]*50/100
    dic.two <- dic.mcmcsamples(res, burnin = b, lik = lik)
    return(dic.two)
}

load("../data/data_for_BiSSE.RData")
load("../data/results_100_phylo_bisse.RData")

state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run the DIC estimate in parallel.
## This is going to run only the first 5 mcmc chains.
## Just to check whether the previous results were ok.
index <- 1:5
tasks <- list(
    job1 <- function() mclapply(index, FUN = function(x) foo.dic.one(mcmc.onerate[[x]][[1]][5000:10000,-1],
																	 tree.genus[[x]], st, unres)
                             , mc.cores = 6),
    job2 <- function() mclapply(index, FUN = function(x) foo.dic.two(mcmc.tworate[[x]][[1]][5000:10000,-1],
																	 tree.genus[[x]], st, unres)
                              , mc.cores = 6)
)

dic.main.out <- mclapply(tasks, function(f) f(), mc.cores = 20)

dic.main.one <- dic.main.out[[1]]
dic.main.two <- dic.main.out[[2]]

save(dic.main.two, dic.main.one, file = "../data/dic_BiSSE_main_results.RData")
