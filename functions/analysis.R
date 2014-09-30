run.bisse <- function(tree, st, unres, steps){
    ## tree = phylo; st = vector of states; unres = unresolved matrix; steps = number of steps of the mcmc chain
    ## save = save every x steps; file = name of the file.
    ## The returning object of this function will be a list with the mcmc run, the lik function and the
    ## prior distribution.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)

    ## This is the analysis with constrained lambda and mu.
    lik.c <- constrain(lik, lambda1~lambda0, mu1~mu0)

    start <- starting.point.bisse(tree)
    prior.c <- make.prior.exponential(1 / 2 * (start[1] - start[3]))

    fit <- find.mle(lik.c, start[argnames(lik.c)])
    tun <- mcmc(lik.c, fit$par, nsteps = 100, w = rep(1,4), lower = 0, print.every = 0, prior = prior.c)
    w <- diff(sapply(tun[2:5], range))
    run.c <- mcmc(lik.c, fit$par, nsteps = steps, w = w, lower = 0, print.every = 0, prior = prior.c)
    return(list(run.c, lik.c, prior.c))
}
