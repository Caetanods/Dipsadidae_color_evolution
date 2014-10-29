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
