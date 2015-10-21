## Prepare data set and phylogenetic trees to run analysis.
## User needs to make sure that the correct data is loaded.

library(diversitree)
library(geiger)

## Load custom functions
source("./functions/data-prepare.R")

## Load trees, species list and color data.
tree <- read.nexus("./data/beast_ingroup_100_posterior_sample.nex")
genus.name <- read.table("./data/species_list.csv", header = TRUE, sep = ",", as.is =TRUE)
full.data <- read.csv("./data/coloration_data.csv", as.is = TRUE)[,-c(1,6)]

## Create the genus level tree.
genera <- genus.name[,-2]
tree.genus <- to.genus.tree(tree, genera)

###################################################
## Load the full data set to make the unresolved table to run BiSSE in the genera tip trees:
## This data set assume that CON (juv) is equal to CON. This means that species in which only
##      the juveniles are aposematic are coded as APO.
## The other polymorphisms are still in this table.

states <- full.data$Category_Std

#####################################################################################################
## DATA FOR MAIN ANALYSIS:

## Transforming the polymorphisms into single categories.
## Change those transformations for alternative codding.

## Aposematic category. Here species with bright colors in the venter are classified as 'contrasting'.
states[which(states == "CON")] <- 1
states[which(states == "VEN")] <- 1
states[which(states == "CON+VEN")] <- 1

## Cryptic category:
states[which(states == "CR")] <- 0
states[which(states == "VIP")] <- 0
states[which(states == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorphic here.
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
one.sp.m <- as.matrix(full.data[full.data[,"Genus"] %in% one.sp,][,c("Genus","BISSE")])
more.sp.m <- cbind(more.sp, NA)
st <- rbind(one.sp.m, more.sp.m)

## Now I need to do the table of diversitree to run BiSSE:
unres <- table(full.data[,c("Genus","BISSE")])
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

#####################################################################################################
## DATA FOR ALTERNATIVE CATEGORIZATION A: VEN -> CR
## In this analysis the species which bright colors are found in the venter only are set as cryptic.

altA <- full.data$Category_Std

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
full.data$altA <- altA

#####################################################################################################
## DATA FOR ALTERNATIVE CATEGORIZATION B: Contrasting juveniles -> CR
## This category is the same as the main analysis. However, all species that were assigned as contrasting
##      based on the juvenile color pattern, now are classified as cryptic.

altB <- full.data$Category_Alt

## Transforming the polymorphisms into single categories.

## Aposematic category:
altB[which(altB == "CON")] <- 1
altB[which(altB == "VEN")] <- 1
altB[which(altB == "CON+VEN")] <- 1

## Changes: contrasting juveniles to the cryptic category (JUV -> CR):
altB[which(altB == "JUV")] <- 0

## Cryptic category:
altB[which(altB == "CR")] <- 0
altB[which(altB == "VIP")] <- 0
altB[which(altB == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorphic here.
altB[which(!altB == "1")] <- 0
full.data$altB <- altB

#####################################################################################################
## DATA FOR ALTERNATIVE CATEGORIZATION C
## This analysis assumes that species which color pattern is similar to coral-snakes are assigned to one
##      category while all the others comprise the other category. Simply put, this analysis is
##      "species that look like coral-snakes" vs. "the rest".

altC <- full.data$Category_Cor

## Fake coral snake category:
altC[which(altC == "C")] <- 1

## Cryptic category. In this case is all the species that are not similar to a coral snake.
altC[which(altC == "N")] <- 0

full.data$altC <- altC

#####################################################################################################
#####################################################################################################

## Prepare the unresolved objects to run the BiSSE MCMC.
## This section of the script will prepare the unresolved objects for all the alternative analysis.
## This follows the same logic performed before with the data for the main analysis.

one.sp.m.altA <- as.matrix(full.data[full.data[,"Genus"] %in% one.sp,][,c("Genus","altA")])
st.altA <- rbind(one.sp.m.altA, more.sp.m)
one.sp.m.altB <- as.matrix(full.data[full.data[,"Genus"] %in% one.sp,][,c("Genus","altB")])
st.altB <- rbind(one.sp.m.altB, more.sp.m)
one.sp.m.altC <- as.matrix(full.data[full.data[,"Genus"] %in% one.sp,][,c("Genus","altC")])
st.altC <- rbind(one.sp.m.altC, more.sp.m)

st.alt <- list(altA = st.altA, altB = st.altB, altC = st.altC)

unresA <- table(full.data[,c("Genus","altA")])
unresA <- cbind(row.names(unresA),unresA)
rownames(unresA) <- NULL
unresB <- table(full.data[,c("Genus","altB")])
unresB <- cbind(row.names(unresB),unresB)
rownames(unresB) <- NULL
unresC <- table(full.data[,c("Genus","altC")])
unresC <- cbind(row.names(unresC),unresC)
rownames(unresC) <- NULL

unres.alt <- list(altA = unresA, altB = unresB, altC = unresC)

## Now we need the total number of species in one of the columns.
for(i in 1:3){
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
save(st.alt, unres.alt, tree.genus, file = "./data/data_for_BiSSE-alt.RData")
