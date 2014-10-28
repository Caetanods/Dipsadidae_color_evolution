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

to.make.rd.state <- function(states, unresolved, freq0, freq1){
	## Function to create new states based on the frequencies of each state.
	## Does not dependend on the phylogeny. So is equivalent to randomizing tips.
    ## states and unresolved objects as in make.bisse.
    ## freq0 and freq1 are frequencies for the 0 and 1 states.
    ## Return list with the random sampled objects. Same structure of the original data.

    ## Change unresolved block:
    rd.unres <- unresolved
    rd.unres$n0 <- NA
    rd.unres$n1 <- NA
    init <- 1

    for(i in 1:length(rd.unres$Nc)){
        nc <- rd.unres$Nc[i]
        rd <- sample(x = c(0,1), size = nc, replace = TRUE, prob = c(st0.freq, st1.freq))
        rd.state <- c(table(rd))
        ifelse(is.na(rd.state[1]), rd.unres$n0[i] <- 0, rd.unres$n0[i] <- as.numeric(rd.state[1]))
        ifelse(is.na(rd.state[2]), rd.unres$n1[i] <- 0, rd.unres$n1[i] <- as.numeric(rd.state[2]))        
        init <- init+nc
    }
    
    ## Change state vector block:
    rd.st <- states
    change <- as.numeric(which(!is.na(states)))
    rd.st[change] <- sample(x = c(0,1), size = length(change), replace = TRUE, prob = c(st0.freq, st1.freq))

    return(list(random.st = rd.st, random.unres = rd.unres))
}
