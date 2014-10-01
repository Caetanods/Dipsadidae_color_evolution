## Prepare dataset and phylogenetic trees to run analysis.
## User needs to make sure that the correct data is loaded.

library(diversitree)
library(geiger)

## Load custom functions
source("./functions/data-prepare.R")

## Load trees, species list and color data.
tree <- read.nexus("./data/beast_ingroup_100_posterior_sample.nex")
genus.name <- read.table("./data/species_list.csv", header = TRUE, sep = ",", as.is =TRUE)
full.data <- read.csv("./data/coloration_data.csv", as.is = TRUE)

## Create the genus level tree.
tree.genus <- to.genus.tree(tree, genus.name)

## Check how the tree look like:
## plot(tree.genus[[1]], direction="upwards"); axisPhylo(side = 4)

###################################################
## Load the full dataset to make the unresolved table to run BiSSE in the genera tip trees:
## This dataset assume that CON (juv) is equal to CON. This means that species in which only
##      the juveniles are aposematic are coded as APO.
## The other polymorphisms are still in this table.

full.data <- full.data[,-c(1,3,5)]

states <- full.data$Category

## Tranforming the polymorphisms into single categories.
## Change those transformations for alternative codding.

## Aposematic category:
states[which(states == "CON")] <- 1
states[which(states == "VEN")] <- 1
states[which(states == "CON+VEN")] <- 1

## Cryptic category:
states[which(states == "CR")] <- 0
states[which(states == "VIP")] <- 0
states[which(states == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorfic here.
## We are codding them as cryptic (this is the conservative interpretation)
states[which(!states == "1")] <- 0
full.data$BISSE <- states

## First we need a vector of states for the genera.
## The genera that will have clade trees need to be assigned as NA.
t.name <- as.matrix(table(full.data$Genus))
t.name <- cbind(rownames(t.name), t.name)
rownames(t.name) <- NULL

one.sp <- t.name[which(t.name[,2] == 1),][,1]
more.sp <- t.name[which(!t.name[,2] == 1),][,1]
one.sp.m <- as.matrix(full.data[full.data[,3] %in% one.sp,][,-1])
more.sp.m <- cbind(more.sp, NA)
st <- rbind(one.sp.m, more.sp.m)

## Now I need to do the table of diversitree to run BiSSE:
unres <- table(full.data[,c(2,3)])
unres <- cbind(row.names(unres),unres)
rownames(unres) <- NULL

## Now we need the total number of species in one of the columns.
mm <- match(unres[,1], t.name[,1])
unres <- cbind(unres, t.name[mm,2])
unres <- unres[which(!unres[,4] == 1),]
colnames(unres) <- c("tip.label","n0","n1","Nc")
unres <- data.frame(unres, stringsAsFactors = FALSE)

## Make sure that the data.frame have the correct format of states:
unres$Nc <- as.numeric(unres$Nc)
unres$n0 <- as.numeric(unres$n0)
unres$n1 <- as.numeric(unres$n1)

## Save the important objects:
save(st, unres, tree.genus, file = "./data/data_for_BiSSE.RData")
