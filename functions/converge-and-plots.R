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

my.medusa <- function (x, time = TRUE, ...){
    ## Function that corrects problem in the plot MEDUSA of geiger.
	## This function plots the MEDUSA results but allow changes in the
	##      parameters of the graphs.
	## Usage is equal to the 'plot.medusa' function in geiger.
    z <- x$zSummary
    shift <- list(cex = 1, bg = "gray", alpha = 0.75, col = "black", 
        lwd = 1)
    phy <- x$cache$phy
    mm <- match(phy$edge[, 2], z[, "dec"])
    edge.color <- z[mm, "partition"]
    plot(phy, edge.color = edge.color, ...)
    if (time) 
        axisPhylo()
    shifts <- z[(idx <- !is.na(z[, "shift"])), "dec"]
    if (length(shifts)) {
        ss <- z[idx, "shift"]
        ww <- character(nrow(phy$edge))
        ww[match(shifts, phy$edge[, 2])] <- ss
        xx <- numeric(nrow(phy$edge))
        xx[phy$edge[, 2] %in% shifts] <- 1
        edgelabels.auteur(NULL, frame = "circle", cex = ifelse(xx == 
            1, shift$cex, 1e-10), pch = ifelse(xx == 1, 21, NA), 
            bg = geiger:::.transparency(shift$bg, shift$alpha), col = shift$col, 
            lwd = shift$lwd)
        edgelabels.auteur(ww, frame = "none", cex = ifelse(xx == 
            1, shift$cex/3, 1e-10), pch = NA, col = shift$col)
    }
}
