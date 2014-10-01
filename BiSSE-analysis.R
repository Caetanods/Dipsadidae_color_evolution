## This is the main script to reproduce the analysis.
## This will run the analysis. Make convergence checks and produce the graphs.

library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run example of analysis:
res.constrain <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis:
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:100
## mcmc.onerate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = st, unres = unres, tun.steps = 100, chain.steps = 100000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = st, unres = unres, tun.steps = 100, chain.steps = 100000, constrain = "FALSE"), mc.cores = 10)

###########################
## Load the results (Comment this part if running the complete analysis)
## Download the file from FigShare.

download.file(url = "http://files.figshare.com/1696849/results_100_phylo_bisse.RData", destfile = ".results_bisse.RData", method = "wget")
load("./results_bisse.RData")

###########################

## Check for convergence using the coda package:

library(coda)

## Create directory to save plots:
dir.create("./mcmc_coda_plots")

## Trace plots
dim(mcmc.onerate[[1]][[1]])[1]
index <- round(seq(1, dim(mcmc.onerate[[1]][[1]])[1], length.out = 500))
runs <- length(mcmc.onerate)

ng <- paste("./mcmc_coda_plots/trace_cons_",seq(1,runs),".pdf", sep="")
for(i in 1:runs){
    pdf(ng[i])
    chain <- mcmc.onerate[[i]][[1]]
    plot(chain$p[index], type = "l")
    dev.off()
}

## Taking out a 25% burn-in:
burn <- list()
for(i in 1:runs){
    burn[[i]] <- mcmc.onerate[[i]][[1]][(dim(res[[i]][[1]])[1]*25/100):dim(res[[i]][[1]])[1],]
}

lb <- paste("./mcmc_coda_plots/lambda_",seq(1,runs),".pdf", sep="")
mu <- paste("./mcmc_coda_plots/mu_",seq(1,runs),".pdf", sep="")
mk <- paste("./mcmc_coda_plots/markov_",seq(1,runs),".pdf", sep="")

for(i in 1:runs){
    pdf(file = lb[i])
    profiles.plot(burn[[i]][c("lambda0")], col.line = c("yellow"))
    dev.off()
    
    pdf(file = mu[i])
    profiles.plot(burn[[i]][c("mu1")], col.line = c("yellow"))
    dev.off()

    pdf(file = mk[i])
    profiles.plot(burn[[i]][c("q01","q10")], col.line = c("yellow","blue"), legend.pos = "topright")
    dev.off()
}

## Now need to run a coda analysis to make sure that this is the converged distribution.

sm <- list()
hei <- list()
for(i in 1:10){
    res <- res.cr[[i]][[1]]
    mc <- res[,c(-1,-8)]
    mc.dat <- as.mcmc(mc)
    sm[[i]] <- summary(mc.dat)
    hei[[i]] <- heidel.diag(mc.dat)
}

hei[[10]]
## Good! We are in stationarity.

## Now the combined posterior with half burn in for all the 10 trees:

burn.comb <- list()
for(i in 1:10){
    burn.comb[[i]] <- res.cr[[i]][[1]][(dim(res.cr[[i]][[1]])[1]*50/100):dim(res.cr[[i]][[1]])[1],]
}
comb <- do.call(rbind, burn.comb)
head(comb)
dim(comb)
comb.one.rate <- comb
save(comb.one.rate, file = "combine.one.rate.RData")

## write.table(comb, file = "Comb.posterior.bisse.cr.txt", sep = ",")

## Now the graphs for the combined posterior:

pdf(file = "comb.lambda.pdf")
profiles.plot(comb[c("lambda0")], col.line = c("blue"))
dev.off()

pdf(file = "comb.mu.pdf")
profiles.plot(comb[c("mu1")], col.line = c("blue"))
dev.off()

pdf(file = "comb.markov.pdf")
profiles.plot(comb[c("q01","q10")], col.line = c("yellow","blue"), legend.pos = "topright")
dev.off()

## pdf(file = "comb.net_diver.pdf")
## profiles.plot(comb[c("lambda0","lambda1")] - comb[c("mu0","mu1")], col.line = c("yellow","blue")
##               , legend.pos = "topright")
## dev.off()

## The results are very similar. Are the trees the same trees?
## First check only the topology:
dt <- matrix(nrow = 10, ncol = 10)
rownames(dt) <- seq(1:10)
colnames(dt) <- seq(1:10)

for(i in 1:10){
    for(j in 1:10){
        dt[i,j] <- dist.topo(tree.gen[[i]],tree.gen[[j]])
    }
}

plot(tree.gen[[3]])
