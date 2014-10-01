to.mcmc.plot <- function(mcmc, run, dir, burn){
    ## mcmc: The log from BiSSE.
    ## run: String to identify the run.
    ## dir: Path for the directory to save plots.
    ## burn: Burnin for the diversitree profile plots.

    dir.create(dir)
    index <- round(seq(1, dim(mcmc)[1], length.out = 500))

    ng <- paste(dir, "trace_", run,".pdf", sep="")
    pdf(ng)
    plot(mcmc$p[index], type = "l")
    dev.off()

    ## Taking out burn-in:
    burn <- mcmc[(dim(mcmc)[1]*25/100):dim(mcmc)[1],]

    lb <- paste(dir, "lambda_", run,".pdf", sep="")
    mu <- paste(dir, "mu_", run,".pdf", sep="")
    mk <- paste(dir, "markov_", run,".pdf", sep="")

    pdf(file = lb)
    profiles.plot(burn[c("lambda0")], col.line = c("yellow"))
    dev.off()
    
    pdf(file = mu)
    profiles.plot(burn[c("mu1")], col.line = c("yellow"))
    dev.off()

    pdf(file = mk)
    profiles.plot(burn[c("q01","q10")], col.line = c("yellow","blue"), legend.pos = "topright")
    dev.off()
}

to.heidel.diag <- function(bisse){
    ## bisse: output of the 'run.bisse' function.
    
    make.heidel <- function(mcmc){
    res <- mcmc[,c(-1,-8)]
    mc.dat <- as.mcmc(res)
    sm <- summary(mc.dat)
    hei <- heidel.diag(mc.dat)
    return(list(summary = sm, diagnostic = hei))
    }

    out <- lapply(1:length(bisse), FUN = function(x) make.heidel(bisse[[x]][[1]]) )
}
