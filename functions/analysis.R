run.bisse <- function(tree, st, unres, tun.steps, chain.steps, constrain = "TRUE"){
    ## tree = phylo; st = vector of states; unres = unresolved matrix; steps = number of steps of the mcmc chain
    ## save = save every x steps; file = name of the file.
    ## The returning object of this function will be a list with the mcmc run, the lik function and the
    ## prior distribution.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
	
	if(constrain == "TRUE"){
		lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
		print("Constrained model.")
	} else { print("Full model.") }

	## Flag:
	print("Start analysis...")

    start <- starting.point.bisse(tree)
    prior <- make.prior.exponential(1 / 2 * (start[1] - start[3]))

    fit <- find.mle(lik, start[argnames(lik)])

	## Flag:
	print("ML estimate finished...")

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = rep(1,4), lower = 0, print.every = 0, prior = prior)

	## Flag:
	print("MCMC tunning finished...")
	print("Starting MCMC chain...")

    w <- diff(sapply(tun[2:5], range))
    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

	## Flag:
	print("Done!")

    return(list(run, lik, prior))
}
