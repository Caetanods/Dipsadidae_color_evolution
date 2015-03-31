## This script will do some simulations to check how the type 1 error affect my results.

## Sim1: Randomizing the data on the tips.
## For this I will create a new dataset of tip data generated with the same frequency of
##     states of the real dataset.

library(multicore)
library(diversitree)
library(phytools)

source("./functions/analysis.R")
load("./data/data_for_BiSSE.RData")

## Create state vector for BiSSE with real data:
state <- as.numeric(st[,2])
names(state) <- st[,1]

## Get total of species and frequency of each state.
## Sum of the 'unres' col plus the monotypic genera.
total.sp <- sum(unres$Nc) + length(which(!is.na(state)))
st1.freq <- ( sum(unres$n1) + length(which(state == "1")) ) / total.sp
st0.freq <- ( sum(unres$n0) + length(which(state == "0")) ) / total.sp

## Function to create the 'unres' and 'state' objects from sampled data with defined frequencies:
rd.data <- to.make.rd.state(state, unres, st0.freq, st1.freq)

## Run example of analysis under the full model.
#res.constrain <- run.bisse(tree = tree.genus[[1]], st = rd.data$random.st, unres = rd.data$random.unres
#                         , tun.steps = 10, chain.steps = 10, constrain = "FALSE")

## Sim2: Simulating a no-effect MK2 model.
## For this we need the complete resolved phylogeny at the species level.

## Create the diver object. Equal to the rich parameter for MEDUSA.
diver <- unres[,c(1,4)]

## Need to create polytomies. This function is too slow!!!
#phy.poly <- lapply(1:length(tree.genus), FUN = function(x) to.make.poly(diver, tree.genus[[x]]) )
#save(phy.poly, file="./data/phy.poly.RData")
load("./data/phy.poly.RData") ## The above line take too much time to run.

## First download and source the Polytomy resolver from the MEE website.
## Ref: Kuhn, T. S., A. Ø. Mooers, and G. H. Thomas. 2011. A simple polytomy resolver for dated phylogenies. Methods in Ecology and Evolution 2:427–436.
download.file(url = "http://onlinelibrary.wiley.com/store/10.1111/j.2041-210X.2011.00103.x/asset/supinfo/MEE3_103_sm_PolytomyResolver.R?v=1&s=8cc376399cf65663d730dd285450339b27c3c1bc", destfile = "./functions/MEE3_103_sm_PolytomyResolver.R", method = "wget")
source("./functions/MEE3_103_sm_PolytomyResolver.R")

## Create directory. Create BEAST files. Run BEAST. Get results:
## Warning: Need BEAST and a Unix (tested in Ubuntu) environment to run:
dir.create(path = "./Beast_files")
dir.create(path = "./Beast_files/xml")
dir.create(path = "./Beast_files/run_results")

sapply(1:length(phy.poly), FUN = function(x) PolytomyResolver(phy.poly[[x]], error = 0.005,file.out = paste0("./Beast_files/xml/phy.mksim.", x) ) )

system(command = "rm ./Beast_files/xml/*.nex")

to.run.BEAST <- function(xml){
    system(
mclapply(1:length(phy.poly),
