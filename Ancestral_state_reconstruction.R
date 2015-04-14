## This is the ancestral state reconstruction integrated over all 100 trees to show that there is no
##         difference in the root state dependent on the phylogenetic uncertainty.

library(diversitree)
library(geiger)

## Get the parameter estimates:
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")
## Get the trees and data:
load("./data/data_for_BiSSE.RData")
## Create state vector for BiSSE:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Mean parameter estimate for each of the 100 trees:
par.one.rate <- t(sapply(1:100, function(x) apply(mcmc.onerate[[x]][[1]][5000:10000,], 2, mean) ))
par.two.rate <- t(sapply(1:100, function(x) apply(mcmc.tworate[[x]][[1]][5000:10000,], 2, mean) ))

## Function to get the MLE ASR with a vector of mean parameter values and a tree:
to.asr <- function(mean.par, phy, st, unres){
    lik <- make.bisse(phy, states = st, unresolved = unres)
    asr <- asr.marginal(lik, mean.par)
    return(asr)
}

## Apply the function to the results:
## Need to fix the one.rate parameter matrix. BiSSE requires the full model.
par.one.rate <- par.one.rate[,c(-1,-6)]
par.one.rate.fix <- par.one.rate[,c(1,1,2,2,3,4)]
colnames(par.one.rate.fix) <- c("lambda0","lambda1","mu0","mu1","q01","q10")

## Run ancestral state reconstruction:
## The output of this function is the nodes as columns and the probability of each state as rows.
## State 0 is in row 1 and state 1 is on row 2.
## asr.one.rate <- lapply(1:100, FUN = function(x) to.asr(par.one.rate.fix[x,], tree.genus[[x]], state, unres))
## asr.two.rate <- lapply(1:100, FUN = function(x) to.asr(par.two.rate[x,c(-1,-8)], tree.genus[[x]], state, unres))
## Save the ancestral state reconstruction for all trees:
## save(asr.one.rate, asr.two.rate, file = "./data/asr_100_phylo.RData")
load("./data/asr_100_phylo.RData")

## Get the distribution of ASR of the root for all the trees:
root.one.rate <- t(sapply(asr.one.rate, FUN = function(x) x[,1]))
colnames(root.one.rate) <- c("State0","State1")
root.two.rate <- t(sapply(asr.two.rate, FUN = function(x) x[,1]))
colnames(root.two.rate) <- c("State0","State1")

## We can see that the root is reconstructed to be in state 1 when the model has trait-dependent
##    diversification and to the state 0 when the diversification is independent of the trait.

#########################################################
## Get the MCC tree:
tree.mcc <- read.tree("./data/ladder.tree.mcc.tre")
tree.mcc <- rescale(tree.mcc, model = "depth", 1.0)
tree.mcc <- multi2di(tree.mcc)

## Make the ancestral state reconstruction under the chosen model. The chosen model was the simpler
##    model where diversification is independent of the traits.
## We are going to use the mean parameter estimate across all trees:
mean.par.one.rate <- apply(par.one.rate.fix, 2, mean)
asr.mcc <- to.asr(mean.par.one.rate, tree.mcc, state, unres)

#########################################################
## Prepare the nice ancestral state reconstruction figure:
mm <- which(!tree.mcc$tip.label %in% unres[,1])
st.bar <- data.frame(tip.label = tree.mcc$tip.label[mm], state = state[tree.mcc$tip.label[mm]])
rownames(st.bar) <- NULL
st.bar <- data.frame(st.bar, n0 = rep(0, times = dim(st.bar)[1]), n1 = rep(0, times = dim(st.bar)[1]))
st.bar$n0[which(st.bar$state == 0)] <- 1
st.bar$n1[which(st.bar$state == 1)] <- 1
st.bar <- st.bar[,-2]

st.bar.total <- unres[,-4]
st.bar.total <- rbind(st.bar.total, st.bar)

## Check if we have all species:
length(tree.mcc$tip.label) == dim(st.bar.total)[1]

mmm <- match(tree.mcc$tip.label, st.bar.total[,1])
tt <- st.bar.total[mmm,2] + st.bar.total[mmm,3]
names(tt) <- tree.mcc$tip.label

st1 <- which(!st.bar.total[mmm,3] == 0 & !st.bar.total[mmm,2] == 0) # Have aposematic and cryptic species.
st11 <- which(st.bar.total[mmm,2] == 0) # Have only aposematic species.
st0 <- which(!st.bar.total[mmm,2] == 0) # Have criptic species.

## Make the pdf for the figure:
pdf("asr_bar_mcc_one.rate.pdf", width = 15, height = 20)
plot.phylo(tree.mcc, font = 4, label.offset = 0.03, cex = 1.0, adj = 1, x.lim = 2, edge.width = 3)
nodelabels(pie = asr.mcc[1,], piecol = c("grey","red"), cex= .35)
axisPhylo(side = 3, pos = c(70,0))

ad.t0 <- 1.275; incr <- 0.007; wd <- 0.2

## Add bars in three sections: All gray, red in top of gray, all red.
for(i in 1:length(st0)){ ## st0
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,2][st0[i]]), ad.t0+(incr*st.bar.total[mmm,2][st0[i]])
                , ad.t0), y = c(st0[i]-wd, st0[i]-wd, st0[i]+wd, st0[i]+wd), col = "grey", lend = 1) }
for(i in 1:length(st1)){ ## st1 on top of st0
    polygon(x = c(ad.t0+(incr*st.bar.total[mmm,2][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]]))
            , y = c(st1[i]-wd, st1[i]-wd, st1[i]+wd, st1[i]+wd), col = "red", lend = 1) }
for(i in 1:length(st11)){ ## st1 only
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,3][st11[i]]), ad.t0+(incr*st.bar.total[mmm,3][st11[i]])
                , ad.t0), y = c(st11[i]-wd, st11[i]-wd, st11[i]+wd, st11[i]+wd), col = "red", lend = 1) }
dev.off()
