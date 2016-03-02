run.bisse <- function(tree, st, unres, tun.steps, chain.steps, constrain = "TRUE", flag = NULL){
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
		print(paste(flag, "- Constrained model.", sep=""))
	} else {
		w.init <- rep(1,6) 
		print(paste(flag,"- Full model.", sep=""))
	}

	## Flag:
	print(paste(flag, "- Start analysis...", sep=""))

    start <- starting.point.bisse(tree)
    prior <- make.prior.exponential(1 / 2 * (start[1] - start[3]))

    fit <- find.mle(lik, start[argnames(lik)])

	## Flag:
	print(paste(flag, "- ML estimate finished...", sep=""))

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = w.init, lower = 0, print.every = 0, prior = prior)

	if(constrain == "TRUE"){
		w <- diff(sapply(tun[2:5], range))
	} else {
		w <- diff(sapply(tun[2:7], range))
	}

	## Flag:
	print(paste(flag, "- Starting MCMC chain...", sep=""))

    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

	write.csv(run, file = paste(flag, "_mcmc.csv", sep=""))
    return(list(run, lik, prior))

	## Flag:
	print(paste(flag, "- Done!", sep=""))
}

run.bisse.mle <- function(tree, st, unres, constrain = "TRUE", flag = NULL){
    ## Function run a maximum estimate of the BiSSE model and return results.
    ## tree = phylo; st = vector of states; unres = unresolved matrix
    ## constrain = use the contrained model; flag = the id of the file and print on screen.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
	
	if(constrain == "TRUE"){
		lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
		print(paste(flag, "- Constrained model.", sep=""))
	} else {
		print(paste(flag,"- Full model.", sep=""))
	}

    start <- starting.point.bisse(tree)

    fit <- find.mle(lik, start[argnames(lik)])

	## Flag:
	print(paste(flag, "- ML estimate finished...", sep=""))

    return(list(likelihood = lik, MLE = fit))
}

run.bisse.bound.mle <- function(tree, st, unres, low, up, constrain = "TRUE", flag = NULL){
    ## Function run a maximum estimate of the BiSSE model and return results.
	## This uses the 'L-BFGS-B' optimization method with bounds to par estimates.
    ## tree = phylo; st = vector of states; unres = unresolved matrix
    ## constrain = use the contrained model; flag = the id of the file and print on screen.
	## low; up = the lower and upper bound for the parameters estimates.
    
    tree <- multi2di(tree)
    mmm <- match(tree$tip.label, names(st))
    lik <- make.bisse(tree, st[mmm], unresolved = unres)
	
	if(constrain == "TRUE"){
		lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
		print(paste(flag, "- Constrained model.", sep=""))
	} else {
		print(paste(flag,"- Full model.", sep=""))
	}

    start <- starting.point.bisse(tree)

    ## fit <- find.mle(lik, start[argnames(lik)], method = "optim"
	## 	   , optim.method = "L-BFGS-B", lower = low, upper = up)
    fit <- find.mle(lik, start[argnames(lik)], method = "optim"
		   , lower = low, upper = up)

	## Flag:
	print(paste(flag, "- ML estimate finished...", sep=""))

    return(list(likelihood = lik, MLE = fit))
}

run.bisse.split <- function(tree, st, unres, nodes, tun.steps, chain.steps, constrain = "TRUE", flag = NULL){
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
	print(paste(flag, "Constrained model.", sep=""))
    } else {
	w.init <- rep(1,12) 
	print(paste(flag, "Full model.", sep=""))
    }

    ## Flag:
    print(paste(flag, "Start analysis...", sep=""))

	## Need also to change the number of parameters for the start point.
    start <- starting.point.bisse(tree)
	start.split <- rep(start, 2)
	names(start.split) <- argnames(lik)
    prior <- make.prior.exponential(1 / 2 * (start.split[1] - start.split[3]))
    fit <- find.mle(lik, start.split[argnames(lik)])

    ## Flag:
    print(paste(flag, "ML estimate finished...", sep=""))

    tun <- mcmc(lik, fit$par, nsteps = tun.steps, w = w.init, lower = 0, print.every = 0, prior = prior)

    if(constrain == "TRUE"){
	w <- diff(sapply(tun[2:9], range))
    } else {
	w <- diff(sapply(tun[2:13], range))
    }

    ## Flag:
    print(paste(flag, "Starting MCMC chain...", sep=""))

    run <- mcmc(lik, fit$par, nsteps = chain.steps, w = w, lower = 0, print.every = 0, prior = prior)

	write.csv(run, file = paste(flag, "_mcmc_split.csv", sep=""))
    return(list(run, lik, prior))

    ## Flag:
    print(paste(flag, "Done!", sep=""))
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
        rd <- sample(x = c(0,1), size = nc, replace = TRUE, prob = c(freq0, freq1))
        rd.state <- c(table(rd))
        ifelse(is.na(rd.state[1]), rd.unres$n0[i] <- 0, rd.unres$n0[i] <- as.numeric(rd.state[1]))
        ifelse(is.na(rd.state[2]), rd.unres$n1[i] <- 0, rd.unres$n1[i] <- as.numeric(rd.state[2]))        
        init <- init+nc
    }
    
    ## Change state vector block:
    rd.st <- states
    change <- as.numeric(which(!is.na(states)))
    rd.st[change] <- sample(x = c(0,1), size = length(change), replace = TRUE, prob = c(freq0, freq1))

    return(list(random.st = rd.st, random.unres = rd.unres))
}

to.make.poly <- function(diver, phy){
    ## Create a phylogeny with polytomies for clade data.
    ## diver is a matrix with first column is the tip.labels and second column is the total
	## Note that this function have custom versions of function in ape and phytools.
	## Customization is to improve speed.
    
    get.faster <- function(phy, tip){
        ## Modification of the ape getMRCA function.
        ## Faster than findMRCA
        tip <- which(phy$tip.label %in% tip)
        Ntip <- length(phy$tip.label)
        seq.nod <- .Call(seq_root2tip, phy$edge, Ntip, phy$Nnode)
        sn <- seq.nod[tip]
        MRCA <- Ntip + 1
        i <- 2
        repeat {
            x <- unique(unlist(lapply(sn, "[", i)))
            if (length(x) != 1) 
                break
            MRCA <- x
            i <- i + 1
        }
        MRCA
    }

    add.species <- function(tree, species, genus = NULL){
        ## This is just a simpler version of 'phytools' function add.species.to.genus.
        ## This function does not have the checkings. As a results is faster.
        ## Original 'phytools' function is more flexible.
    
        ## The 'where' parameter here is considered always == 'root'.
        ii <- grep(paste(genus, "_", sep = ""), tree$tip.label)
        if (length(ii) > 1) {
            nn <- get.faster(tree, tree$tip.label[ii])
            tree <- bind.tip(tree, gsub(" ", "_", species), where = nn)
        } else {
            nn <- ii
            tree <- bind.tip(tree, gsub(" ", "_", species), where = nn, 
                             position = 0.5 * tree$edge.length[which(tree$edge[, 2] == nn)])
        }
        return(tree)
    }
    
    phy$tip.label <- paste(phy$tip.label,"_sp1",sep="")
    sp.list <- lapply(1:dim(diver)[1], FUN = function(x) paste(diver[x,1],"_sp",2:diver[x,2],sep=""))
    for(i in 1:length(sp.list)){
        print(paste("Clade ",i," of ",length(sp.list)," with ",diver[i,2]," species.",sep=""))
        for(j in 1:length(sp.list[[i]])){
            phy <- add.species(phy, species=sp.list[[i]][j], genus=diver[i,1])
        }
    }
    return(phy)
}

dic.mcmcsamples <- function(x, burnin=0, lik){
  ## Compute deviance information criterion from mcmcsamples
  ## This is the version corrected. It will work with the RCran version
  ##  of 'diversitree'.

  ## if (!inherits(x, "mcmcsamples"))
  ##  stop("this function is only designed for diversitree's mcmcsamples object")
  ## lik <- attr(x, "func")
  ## p <- coef(x, burnin=burnin)
  p <- x  ## Because I am reading from a csv. So 'x' is a matrix with parameter values.
  dev <- -2 * apply(p, 1, lik)
  ## The previous line seems not to be needed. Note that the posterior produced by the BiSSE
  ##   analysis already have the likelihood of each generation of the mcmc. This means
  ##   that this line is nothing but a recalculation of this quantity.

  ## estimate effective number of parameters
  dbar <- mean(dev)

  ## deviance of posterior means:
  post.means <- colMeans(p)

  ## evaluate deviance at the mean posterior estimate
  dhat <- -2 * lik(post.means)
  pd <- dbar - dhat

  ## calculate dic
  dic <- dbar + pd
  unname(dic)
}

dic.lite <- function(x, p, burnin=0, lik){
  ## Compute deviance information criterion from mcmcsamples.
  ## This version does not recalculates the likelihood for entire joined posterior distribution.
  ## x = matrix with only the parameters of the posterior. Other columns need to be excluded.
  ## p = The log likelihood as present in the MCMC results from the BiSSE analyses.
  ## burnin = optional.
  ## lik = The likelihood function to calculate the likelihood of the mean parameter values.

  ## Get the deviance from the posterior
  dev <- -2 * p

  ## estimate effective number of parameters
  dbar <- mean(dev)

  ## deviance of posterior means:
  post.means <- colMeans(x)

  ## evaluate deviance at the mean posterior estimate
  dhat <- -2 * lik(post.means)
  pd <- dbar - dhat

  ## calculate dic
  dic <- dbar + pd
  unname(dic)
}

sims <- function(pars, time, freq.root){
    ## Function always return a tree with given parameters.
    ## This function will stop the BiSSE tree simulation using time only.
    ## The root state is drawn for a binomial distribution and parameters
    ##     from a given distribution.
    ## Arguments:
    ## pars: data frame with 6 columns lambda0; lambda1; mu0; mu1; q01; q10.
    ## time: numeric. the depth of the simulated tree.
    ## freq.root: vector. frequency of the state 1 and 0 in the root.
    repeat
        {
            ll <- sample(1:dim(pars)[1], 1)
            par <- as.numeric(pars[ll,])
            root <- sample(c(0,1), size = 1, prob = freq.root)
            phy <- tree.bisse(par, max.t = time, x0 = root)
            if(!is.null(phy))
                {
                    break
                }
        }
    return(phy)
}

get.medusa.rates <- function(medusa.res, phy){
    ## Function to summarize the results from the MEDUSA analysis.
    ## This function produces a table with the tip names and the rate associated with
    ## each tip.
    ## medusa.res = res object from medusa
    ## phy = phylogenetic tree

    model <- cbind(medusa.res[[3]][[4]], medusa.res[[3]][[1]])
    
    tt.l <- list()
    ss <- vector()

    for(i in 1:dim(model)[1] ){
        tt <- tips(tree.genus[[2]], as.numeric(model[i,1]) )
        ss[i] <- length(tt)
        tt.m <- cbind(tt, model[i,2])
        tt.l[[i]] <- tt.m
        res <- list(shifts = tt.l, length = ss)
    }

    ord <- order(res$length)
    for( j in 1:(length(ord)-1) ){
        for( k in (j+1):length(ord) ){
            index <- res[[ 1 ]][[ ord[k] ]][ , 1 ] %in% res[[ 1 ]][[ ord[j] ]][ , 1 ]
            if( sum(index) != 0 ){
                res[[ 1 ]][[ ord[k] ]][index , 2] <- res[[ 1 ]][[ ord[j] ]][ , 2]
            } else { }
        }
    }
    
    res.t <- res[[ 1 ]][[ ord[length(ord)] ]]
    out <- as.numeric(res.t[,2])
    names(out) <- res.t[,1]    

    return( out )
}
