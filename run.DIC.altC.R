require(diversitree)
require(parallel)

source("./functions/analysis.R")

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

load("./data/data_for_BiSSE-alt.RData")
load("./data/results_bisse_mcmc_alt_ABC.RData")

stateC <- as.numeric(st.alt[[3]][,2])
names(stateC) <- st.alt[[3]][,1]

## Run the DIC estimate in parallel.
index <- 1:5
tasks <- list(
    job1 <- function() mclapply(index, FUN = function(x) foo.dic.one(mcmc.onerate.C[[x]][5000:10000,c(3:6)],
																	 tree.genus[[x]], stateC,
																	 unres.alt[[3]])
                             , mc.cores = 6),
    job2 <- function() mclapply(index, FUN = function(x) foo.dic.two(mcmc.tworate.C[[x]][5000:10000,c(3:8)], 
																	 tree.genus[[x]], stateC, 
																     unres.alt[[3]])
                              , mc.cores = 6)
)

dic.altC.out <- mclapply(tasks, function(f) f(), mc.cores = 20)

dic.altC.one <- dic.altC.out[[1]]
dic.altC.two <- dic.altC.out[[2]]

save(dic.altC.two, dic.altC.one, file = "./data/dic_BiSSE_altC_results.RData")
