## Prepare dataset and phylogenetic trees to run analysis.
## User needs to make sure that the correct data is loaded.

## Directory of caetano@129.101.170.113 machine:
## ~/Documents/Snakes/diversitree/.
## This machine is using the 3.0.0 version of R.
## This script is using trees sampled from the posterior distribution of the beast analysis.

require(diversitree)

## Load trees and species list
tree <- read.nexus("./data/beast_ingroup_100_posterior_sample.nex")
genus.name <- read.table("./data/species_list.csv", header = TRUE, sep = ",", as.is =TRUE)

## Create the genus level tree.

to.genus.tree <- function(phy, genera){
    ## Picking one species per genus randomly:
    ## phy = multiphylo object
    ## genera = vector with genera names
    genus.list <- list()
    temp <- list()
    gu <- unique(genera)
    for(i in 1:length(gu)){
        temp[[i]] <- genera[which(genera == gu[i])]
        a <- sample(1:dim(temp[[i]])[1], 1)
        genus.list[[i]] <- temp[[i]][a,]
    }
    genus.list <- do.call(rbind,genus.list)

## Need to remember to check if the data table has the same list of genera. If not, the genera need to be
## pruned of the data table. Since it is not possible to place the genus in the phylogeny in the correct place.
genus.list

## Prune the tree only to keep those species:
to.keep <- genus.list[,3]
to.prune <- tree[[1]]$tip.label[!tree[[1]]$tip.label %in% to.keep]
length(to.keep) + length(to.prune) == length(tree[[1]]$tip.label) #Test to see if we got the right tips.
tree.genus <- lapply(1:100, FUN = function(x) drop.tip(tree[[x]], tip=to.prune))

## Note that the depth of the tree is not scale to 1:
plot(tree.genus[[1]], direction="upwards"); axisPhylo(side = 4)

## Checking if all the trees have the same distance from root to tip:
tb <- vector()
for(i in 1:100) tb[i] <- vcv.phylo(tree.genus[[i]])[1,1]
tb
## Note that the trees have almost the same length from root to tip. But not the same.
require(geiger)
tree.genus <- lapply(tree.genus, FUN = function(x) rescale(x, model = "depth", 1.0))
plot.phylo(tree.genus[[1]],direction = "upwards"); axisPhylo(side = 4)

## Change the names of the species to the name of the genera:
tree.gen <- tree.genus
mm <- lapply(1:100, FUN = function(x) match(tree.gen[[x]]$tip.label, genus.list[,3]))
for(i in 1:100){
    tree.gen[[i]]$tip.label <- genus.list[mm[[i]],1]
}
plot(tree.gen[[1]], direction="upwards"); axisPhylo(side = 4)

###################################################
## Use the full dataset to make the unresolved table to run BiSSE in the genera tip trees:
## This dataset assume that CON (juv) is equal to CON. The other polymorphisms are still here.
full.data <- read.csv("dados_coloracao_paulo_full_Jan_23.csv", as.is = TRUE)
full.data <- full.data[,-c(1,4)]

## ####################################################################
## I am considering here the analysis of CON+VEN vs. CR+VIP
## For the analysis of CON vs. VEN+CR+VIP see the last version of the analysis the version 6.

## Check the states of the table:
emp.states <- full.data$Category
unique(emp.states)

emp.states[which(emp.states == "CON")] <- 1
emp.states[which(emp.states == "VEN")] <- 1
emp.states[which(emp.states == "CON+VEN")] <- 1

## Those are for sure "CR":
emp.states[which(emp.states == "CR")] <- 0
emp.states[which(emp.states == "VIP")] <- 0
emp.states[which(emp.states == "CR+VIP")] <- 0

## Now I am going to code for two scenarios.
## First all polymorphisms are going to be "CON+VEN":
states.ven <- emp.states
states.ven[which(!states.ven == 0)] <- 1
full.data$BISSE.ven <- states.ven
## Then all polymorphisms are goin to be "CR":
states.cr <- emp.states
states.cr[which(!states.cr == 1)] <- 0
full.data$BISSE.cr <- states.cr

## First I need a vector of states for the genera.
## The genera that will have clade trees need to be assigned as NA.
t.name <- as.matrix(table(full.data$Genus))
t.name <- cbind(rownames(t.name), t.name)
rownames(t.name) <- NULL

one.sp <- t.name[which(t.name[,2] == 1),][,1]
more.sp <- t.name[which(!t.name[,2] == 1),][,1]
one.sp.m <- as.matrix(full.data[full.data[,2] %in% one.sp,][,-1])
more.sp.m <- cbind(more.sp, NA, NA)
st <- rbind(one.sp.m, more.sp.m)

## Now I need two vectors to use as entry in BiSSE. One to the all "CON" and the other to all "CR".
st.ven <- as.numeric(st[,2])
names(st.ven) <- st[,1]
st.cr <- as.numeric(st[,3])
names(st.cr) <- st[,1]

## Now I need to do the table of diversitree to run BiSSE:
## One for the all "CON" and the other for the all "CR":
unres.ven <- table(full.data[,c(2,3)])
unres.ven <- cbind(row.names(unres.ven),unres.ven)
rownames(unres.ven) <- NULL

unres.cr <- table(full.data[,c(2,4)])
unres.cr <- cbind(row.names(unres.cr),unres.cr)
rownames(unres.cr) <- NULL

## Cannot sum the rows because there are some NA values. Need to get the previous table,
## with the number of species in each clade.
mm <- match(unres.cr[,1], t.name[,1])
unres.cr <- cbind(unres.cr, t.name[mm,2])
unres.cr <- unres.cr[which(!unres.cr[,4] == 1),]

mm <- match(unres.ven[,1], t.name[,1])
unres.ven <- cbind(unres.ven, t.name[mm,2])
unres.ven <- unres.ven[which(!unres.ven[,4] == 1),]

colnames(unres.cr) <- c("tip.label","n0","n1","Nc")
unres.cr <- data.frame(unres.cr, stringsAsFactors = FALSE)
colnames(unres.ven) <- c("tip.label","n0","n1","Nc")
unres.ven <- data.frame(unres.ven, stringsAsFactors = FALSE)

## Make sure that the data.frame have the correct format of states:
unres.cr$Nc <- as.numeric(unres.cr$Nc)
unres.cr$n0 <- as.numeric(unres.cr$n0)
unres.cr$n1 <- as.numeric(unres.cr$n1)

unres.ven$Nc <- as.numeric(unres.ven$Nc)
unres.ven$n0 <- as.numeric(unres.ven$n0)
unres.ven$n1 <- as.numeric(unres.ven$n1)

save.image("beast.mcmc.data_v8_90.RData")
