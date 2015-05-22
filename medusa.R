## This analysis use MEDUSA to search for nodes where net diversification rates change.

## library(multicore)
## 'multicore' will not work anymore. Need to work with parallel::mclapply()
## So the parallel package is already loaded in the new version of R.
library(geiger)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Make richness table for MEDUSA.
total <- unres[,c(1,4)]
single <- tree.genus[[1]]$tip.label[which(!tree.genus[[1]]$tip.label %in% total$tip.label)]
single <- cbind(single,1)
colnames(single) <- colnames(total) <- c("taxon","n.taxa")
rich <- rbind(total, single)

## Constraining the MEDUSA model to find at max two partitions in the phylogenies.
## shifts <- mclapply(1:100, FUN = function(x) medusa(tree.genus[[x]], richness = rich, criterion = "aicc")
##                  , mc.cores = 20)

## save(shifts, file = "medusa.RData")

######################################
## Get results:
load("./data/medusa.RData")

## Scratch for function to summarize result from medusa across phylos.

get.medusa.rates <- function(medusa.res, phy){
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

## Now check the Medusa results across all trees:

## Filter the medusa results to take out the ones that did not converged:
## Index has only analysis which worked fine.
not <- which(sapply(shifts, length) == 1)
index <- 1:100
index <- index[-not]

res <- lapply(index, function(x) get.medusa.rates(shifts[[x]], tree.genus[[x]]) )
res <- do.call(cbind, res)
dim(res)

md <- apply(res, 1, median)
mn <- apply(res, 1, min)
mx <- apply(res, 1, max)
nm <- rownames(res)

pdf("Summarize_medusa.pdf", width = 14)
plot(x=1:69, y = md, col = "blue", ylim = c(0,50), pch = 4
   , axes = FALSE, ylab = "Medusa - yule rate", xlab = "")
axis(side = 2)
text(1:69, par("usr")[3]-0.25, srt = 60, adj= 1, xpd = TRUE, labels = nm, cex=0.8)
points(x=1:69, y = mn, pch = 20)
points(x=1:69, y = mx, pch = 20)
for(i in 1:69) lines(x = c(i,i), y = c(mn[i],mx[i]), lty = 2 )
dev.off()

pdf("Summarize medusa boxblot.pdf", width = 14)
boxplot(t(res), axes = FALSE, ylab = "Medusa - yule rate", xlab = "")
axis(side = 2)
text(1:69, par("usr")[3]-0.25, srt = 60, adj= 1, xpd = TRUE, labels = nm, cex=0.8)
dev.off()

## Plotting medusa results
source("medusa.plot.R")
my.medusa(shifts[[5]], edge.width = 2, )
## This is a corrected version of medusa. Need to push to geiger.
