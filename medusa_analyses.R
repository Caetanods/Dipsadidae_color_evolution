## This analysis use MEDUSA to search for nodes where net diversification rates change.
## Medusa is agnostic to the traits under study and only uses the information of the
##        topology, branch lengths and diversity.

library(parallel)
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

## This runs the MEDUSA analysis with the true diversity of the group and using the
##      aicc criterion to evaluate the jumps in rate under a BD model.
## WARNING: This analysis takes time to run and uses multiple cores.
## Uncomment the following block of code to run the analysis or skip to load the results.

## shifts <- mclapply(1:100, FUN = function(x) medusa(tree.genus[[x]], richness = rich, criterion = "aicc")
##                  , mc.cores = 20)
## save(shifts, file = "medusa.RData")

######################################
## Get results:
load("./data/medusa.RData")

## Now check the Medusa results across all trees:

## Filter the medusa results to take out the ones that did not converged:
## Index has only analysis which worked fine.
not <- which(sapply(shifts, length) == 1)
index <- 1:100
index <- index[-not]
res <- lapply(index, function(x) get.medusa.rates(shifts[[x]], tree.genus[[x]]) )
res <- do.call(cbind, res)

md <- apply(res, 1, median)
mn <- apply(res, 1, min)
mx <- apply(res, 1, max)
nm <- rownames(res)

## Get same order from the MCC in the Figure 1 of manuscript:
## Load the ladderized MCC tree:
mcc.tr <- read.tree("./data/ladder.tree.mcc.tre")
spp <- mcc.tr$tip.label
group <- rep("grey", times = length(spp) )
green <- c("Uromacer","Arrhyton","Haitiophis","Magliophis","Alsophis","Borikenophis"
          ,"Schwartzophis","Hypsirhynchus","Antillophis","Cubophis","Caraiba","Darlingtonia"
          ,"Ialtris") ## 2ca02cff
blue <- c("Hydrops","Pseudoeryx") ## 0055d4ff
red <- c("Pseudotomodon","Tachymenis","Tomodon","Thamnodynastes","Ptychophis") ## d40000ff

group[which(spp %in% green)] <- "#2ca02cff"
group[which(spp %in% blue)] <- "#0055d4ff"
group[which(spp %in% red)] <- "#d40000ff"

## Margin: "bottom, left, top and right"
par( mar=c(5.0, 5.0, 0.5, 1.0) )
boxplot(t(res[spp,]), axes = FALSE, ylab = "", xlab = "Medusa - yule rate"
      , outpch = 20, outcex = 0.8, horizontal = "TRUE", col = group)
axis(side = 1)
## text(1:69, par("usr")[3]-0.25, srt = 60, adj= 1, xpd = TRUE, labels = rownames(res[spp,]), cex=0.8)
text(par("usr")[3]+1.5, 1:69, adj= 1, xpd = TRUE, labels = rownames(res[spp,]), cex=0.8)

## Plotting medusa results on the MCC tree.
## 'my.medusa' is a function which corrects bugs from the original plotting function in 'geiger'.
my.medusa( shifts[[5]], edge.width = 2 )
