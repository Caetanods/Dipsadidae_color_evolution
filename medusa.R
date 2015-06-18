## This analysis use MEDUSA to search for nodes where net diversification rates change.

## library(multicore)
## 'multicore' will not work anymore. Need to work with parallel::mclapply()
## So the parallel package is already loaded in the new version of R.
library(geiger)

source("./functions/analysis.R")
source("./functions/converge-and-plots.R")
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

## Get same order from the MCC in the Figure 1 of manuscript:
## Load the ladderized tree:
mcc.tr <- read.tree("./data/ladder.tree.mcc.tre")

pdf("mcc.tree.pdf", width = 7, height = 14)
plot.phylo(mcc.tr, no.margin = TRUE)
dev.off()

spp <- mcc.tr$tip.label ## Need to get the reverse of the string vector.

## Create group of colors:

group <- rep("grey", times = length(spp) )

green <- c("Uromacer","Arrhyton","Haitiophis","Magliophis","Alsophis","Borikenophis"
          ,"Schwartzophis","Hypsirhynchus","Antillophis","Cubophis","Caraiba","Darlingtonia"
          ,"Ialtris") ## 2ca02cff
blue <- c("Hydrops","Pseudoeryx") ## 0055d4ff
red <- c("Pseudotomodon","Tachymenis","Tomodon","Thamnodynastes","Ptychophis") ## d40000ff

group[which(spp %in% green)] <- "#2ca02cff"
group[which(spp %in% blue)] <- "#0055d4ff"
group[which(spp %in% red)] <- "#d40000ff"

pdf("Summarize medusa boxblot.pdf", width = 7, height = 14)
## Margin: "bottom, left, top and right"
par( mar=c(5.0, 5.0, 0.5, 1.0) )
boxplot(t(res[spp,]), axes = FALSE, ylab = "", xlab = "Medusa - yule rate"
      , outpch = 20, outcex = 0.8, horizontal = "TRUE", col = group)
axis(side = 1)
## text(1:69, par("usr")[3]-0.25, srt = 60, adj= 1, xpd = TRUE, labels = rownames(res[spp,]), cex=0.8)
text(par("usr")[3]+1.5, 1:69, adj= 1, xpd = TRUE, labels = rownames(res[spp,]), cex=0.8)
dev.off()

## Plotting medusa results
my.medusa(shifts[[5]], edge.width = 2 )
## This is a corrected version of medusa. Need to push to geiger.
