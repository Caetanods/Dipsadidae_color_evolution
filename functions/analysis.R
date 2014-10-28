run.bisse <- function(tree, st, unres, tun.steps, chain.steps, constrain = "TRUE"){
    ## tree = phylo; st = vector of states; unres = unresolved matrix; steps = number of steps of the mcmc chain
    ## save = save every x steps; file = name of the file.
    ## The returning object of this function will be a list with the mcmc run, the lik function and the
    ## prior distribution.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
	
	if(constrain == "TRUE"){
		w.init <- rep(1,4)
		lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
		print("Constrained model.")
	} else {
		w.init <- rep(1,6) 
		print("Full model.")
	}

	## Flag:
	print("Start analysis...")

    start <- starting.point.bisse(tree)
    prior <- make.prior.exponential(1 / 2 * (start[1] - start[3]))

    fit <- find.mle(lik, start[argnames(lik)])

	## Flag:
	print("ML estimate finished...")

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = w.init, lower = 0, print.every = 0, prior = prior)

	if(constrain == "TRUE"){
		w <- diff(sapply(tun[2:5], range))
	} else {
		w <- diff(sapply(tun[2:7], range))
	}

	## Flag:
	print("MCMC tunning finished...")
	print("Starting MCMC chain...")

    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

	## Flag:
	print("Done!")

    return(list(run, lik, prior))
}

run.bisse.split <- function(tree, st, unres, nodes, tun.steps, chain.steps, constrain = "TRUE"){
    ## tree = phylo; st = vector of states; unres = unresolved matrix; steps = number of steps of the mcmc chain
    ## save = save every x steps; file = name of the file; node = nodes where to split the model.
    ## The returning object of this function will be a list with the mcmc run, the lik function and the
    ## prior distribution.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))

    ## Prepare the split model:
    lik <- make.bisse.split(tree, st[mmm], unresolved = unres, nodes = nodes, split.t = "Inf")

	## Note that the number of parameters here double at each split of the model.
	## Each run of MEDUSA found exactly one split point. -> 12 parameters total.
    if(constrain == "TRUE"){
	## w.init = free model minus the number of constrained paramters.
	w.init <- rep(1,8)
	lik <- constrain(lik, lambda1.1~lambda0.1, mu1.1~mu0.1, lambda1.2~lambda0.2, mu1.2~mu0.2)
	print("Constrained model.")
    } else {
	w.init <- rep(1,12) 
	print("Full model.")
    }

    ## Flag:
    print("Start analysis...")

	## Need also to change the number of parameters for the start point.
    start <- starting.point.bisse(tree)
	start.split <- rep(start, 2)
	names(start.split) <- argnames(lik)
    prior <- make.prior.exponential(1 / 2 * (start.split[1] - start.split[3]))
    fit <- find.mle(lik, start.split[argnames(lik)])

    ## Flag:
    print("ML estimate finished...")

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = w.init, lower = 0, print.every = 0, prior = prior)

    if(constrain == "TRUE"){
	w <- diff(sapply(tun[2:9], range))
    } else {
	w <- diff(sapply(tun[2:13], range))
    }

    ## Flag:
    print("MCMC tunning finished...")
    print("Starting MCMC chain...")

    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

    ## Flag:
    print("Done!")

    return(list(run, lik, prior))
}
