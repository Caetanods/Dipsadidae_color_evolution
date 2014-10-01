## This script will run the alternative categorization for the colors of the Dipsadidae snakes.
## There are two alternative categorizations:
## a) Species which show ventral contrasting colors as "CR".
## b) Species which only the juvenile show contrastiong colors as "CR".

library(diversitree)

## Load the data and the functions:
source("./functions/data-prepare.R")
source("./functions/analysis.R")
data <- read.csv("./data/coloration_data.csv", as.is = TRUE)[,-c(1,5)]
load("./data/data_for_BiSSE.RData")
rm(st, unres)

###########################
## Alternative categorization A: VEN -> CR

altA <- data$Category_Std

## Tranforming the polymorphisms into single categories.
## Change those transformations for alternative codding.

## Aposematic category:
altA[which(altA == "CON")] <- 1

## Changes: VEN -> CR:
altA[which(altA == "VEN")] <- 0
altA[which(altA == "CON+VEN")] <- 0

## Cryptic category:
altA[which(altA == "CR")] <- 0
altA[which(altA == "VIP")] <- 0
altA[which(altA == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorfic here.
altA[which(!altA == "1")] <- 0
data$altA <- altA
###########################

###########################
## Alternative categorization B: juvenile only -> CR

altB <- data$Category_Alt

## Transforming the polymorphisms into single categories.

## Aposematic category:
altB[which(altB == "CON")] <- 1
altB[which(altB == "VEN")] <- 1
altB[which(altB == "CON+VEN")] <- 1

## Changes: JUV -> CR:
altA[which(altA == "JUV")] <- 0

## Cryptic category:
altB[which(altB == "CR")] <- 0
altB[which(altB == "VIP")] <- 0
altB[which(altB == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorfic here.
altB[which(!altB == "1")] <- 0
data$altB <- altB

###########################

###########################
## Prepare the unresolved data table for BiSSE.
## Do it for the two alternative coddings.

## First we need a vector of states for the genera.
## The genera that will have clade trees need to be assigned as NA.
t.name <- as.matrix(table(data$Genus))
t.name <- cbind(rownames(t.name), t.name)
rownames(t.name) <- NULL

one.sp <- t.name[which(t.name[,2] == 1),][,1]
more.sp <- t.name[which(!t.name[,2] == 1),][,1]
one.sp.m <- as.matrix(data[data[,3] %in% one.sp,][,-c(1,2)])
more.sp.m <- cbind(more.sp, NA, NA)
st <- rbind(one.sp.m, more.sp.m)

## Now I need to do the table of diversitree to run BiSSE.
## One for each category alternative.
unresA <- table(data[,c(3,4)])
unresA <- cbind(row.names(unresA),unresA)
rownames(unresA) <- NULL

unresB <- table(data[,c(3,5)])
unresB <- cbind(row.names(unresB),unresB)
rownames(unresB) <- NULL

unres.alt <- list(altA = unresA, altB = unresB)

## Now we need the total number of species in one of the columns.
for(i in 1:2){
    mm <- match(unres.alt[[i]][,1], t.name[,1])
    unres.alt[[i]] <- cbind(unres.alt[[i]], t.name[mm,2])
    unres.alt[[i]] <- unres.alt[[i]][which(!unres.alt[[i]][,4] == 1),]
    colnames(unres.alt[[i]]) <- c("tip.label","n0","n1","Nc")
    unres.alt[[i]] <- data.frame(unres.alt[[i]], stringsAsFactors = FALSE)
    unres.alt[[i]]$Nc <- as.numeric(unres.alt[[i]]$Nc)
    unres.alt[[i]]$n0 <- as.numeric(unres.alt[[i]]$n0)
    unres.alt[[i]]$n1 <- as.numeric(unres.alt[[i]]$n1)
}

## Save the important objects:
save(st, unres.alt, tree.genus, file = "./data/data_for_BiSSE-alt.RData")

###########################
## Run the analysis:
## First alternative category A:

stateA <- as.numeric(st[,2])
names(stateA) <- st[,1]

## Run example of analysis:
res.constrain.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.A <- run.bisse(tree = tree.genus[[1]], st = stateA, unres = unres.alt[[1]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.A <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateA, unres = unres.alt[[1]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

## Now category B:
stateB <- as.numeric(st[,3])
names(stateB) <- st[,1]

## Run example of analysis:
res.constrain.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "TRUE")
res.free.B <- run.bisse(tree = tree.genus[[1]], st = stateB, unres = unres.alt[[2]], tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Run the complete analysis (with only the first 5 trees):
## WARNING: It depends on 'multicore' and takes time to run.
## library(multicore)
## index <- 1:5
## mcmc.onerate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "TRUE"), mc.cores = 10)
## mcmc.tworate.B <- mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateB, unres = unres.alt[[2]], tun.steps = 100, chain.steps = 10000, constrain = "FALSE"), mc.cores = 10)

