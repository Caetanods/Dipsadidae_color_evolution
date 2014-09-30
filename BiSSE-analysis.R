## This is the main script to reproduce the analysis.
## This will run the analysis. Make convergence checks and produce the graphs.

library(diversitree)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Run example of analysis:
res <- run.bisse(tree = tree.genus[[1]], st = state, unres = unres, steps = 100)

## Run the complete analysis:
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:100
## res <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = st, unres = unres, steps = 10000), mc.cores = 10)

###########################

## Check for convergence using the coda package:

library(coda)

## Create directory to save plots:
dir.create("./mcmc_coda_plots")

## Trace plots
dim(res.cr[[1]][[1]])[1]
index <- round(seq(1, dim(res.cr[[1]][[1]])[1], length.out = 500))
runs <- length(res)

ng <- paste("Trace_",seq(1,runs),".pdf", sep="")
for(i in 1:runs){
    pdf(ng[i])
    chain <- res[[i]][[1]]
    plot(chain$p[index], type = "l")
    dev.off()
}

## Taking out a 25% burn-in:
burn <- list()
for(i in 1:runs){
    burn[[i]] <- res[[i]][[1]][(dim(res[[i]][[1]])[1]*25/100):dim(res[[i]][[1]])[1],]
}

lb <- paste0("lambda_",seq(1,runs),".pdf")
mu <- paste0("mu_",seq(1,runs),".pdf")
mk <- paste0("markov_",seq(1,runs),".pdf")
nt <- paste0("net_diver_",seq(1,runs),".pdf")

for(i in 1:10){
    pdf(file = lb[i])
    profiles.plot(burn[[i]][c("lambda0")], col.line = c("yellow"))
    dev.off()
    
    pdf(file = mu[i])
    profiles.plot(burn[[i]][c("mu1")], col.line = c("yellow"))
    dev.off()

    pdf(file = mk[i])
    profiles.plot(burn[[i]][c("q01","q10")], col.line = c("yellow","blue"), legend.pos = "topright")
    dev.off()

    ## pdf(file = "turn_over.pdf")
    ## profiles.plot(burn[c("mu0","mu1")] / burn[c("lambda0","lambda1")], col.line = c("yellow","blue")
    ##               , legend.pos = "topright")
    ## dev.off()

    ## pdf(file = nt[i])
    ## profiles.plot(burn[[i]][c("lambda0","lambda1")] - burn[[i]][c("mu0","mu1")], col.line = c("yellow","blue")
    ##               , legend.pos = "topright")
    ## dev.off()
}

## Note that the markov transition rate q10 is 10 times higher than the contrary.
## We would expect then that in the long run most of the species would be in the state 0.
## length(states.cr)
## unique(states.cr)
## sum(as.numeric(states.cr), na.rm = TRUE)
## length(states.cr) - sum(as.numeric(states.cr), na.rm = TRUE)

## ## This is a way to test the hypothesis here. We can set a new distribution of the differences off
## ## the net diversification rates dependent of the two traits and see where the expected value will fall.
## net <- burn[c("lambda0","lambda1")] - burn[c("mu0","mu1")]
## colnames(net) <- c("net.diver0","net.diver1")

## ## Compute the cumulative distribution function:
## cdf <- ecdf(net[,2] - net[,1])

## pdf("net_diver_cdf.pdf")
## plot(cdf, main = "CDF of the net diversification of CON - CR")
## dev.off()

## ## This means that the probability of getting values smaller or equal to 0 is smaller than 0.05.
## ## We reject the hypothesis of no difference between the net diversification values with CI equal to 0.95.
## quantile(cdf, probs = 0.05)

## Now need to run a coda analysis to make sure that this is the converged distribution.
require(coda)

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
