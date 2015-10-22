## This script is going to run the alternative category 3.
## The analysis is running in 'macgyver'. ~/Documents/snake_manuscript/snake.colors.of.deception

library(diversitree)
library(parallel)
library(methods)
source("./functions/analysis.R")
load("./data/data_for_BiSSE-alt.RData")

stateC <- as.numeric(st.alt[[3]][,2])
names(stateC) <- st.alt[[3]][,1]

## Run the complete analysis (with only the first 5 trees):
index <- 1:5

tasks <- list(
    job1 <- function() mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateC
                                                               , unres = unres.alt[[3]], tun.steps = 100
                                                               , chain.steps = 10000, constrain = "TRUE"
                                                               , flag = paste("stB_constrain_",x,sep=""))
                             , mc.cores = 6),
    job2 <- function() mclapply(index, FUN = function(x) run.bisse(tree = tree.genus[[x]], st = stateC
                                                               , unres = unres.alt[[3]], tun.steps = 100
                                                               , chain.steps = 10000, constrain = "FALSE"
                                                               , flag = paste("stB_free_",x,sep=""))
                              , mc.cores = 6)
)

altC.out <- mclapply(tasks, function(f) f(), mc.cores = 20)

save(altC.out, file = "altC_out.RData")
