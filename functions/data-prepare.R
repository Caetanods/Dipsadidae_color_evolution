to.genus.tree <- function(phy, genera){
    ## Picking one species per genus randomly:
    ## phy = multiphylo object
    ## genera = data matrix. first column with the genera names. second column with species names. This is the same species of the original phylogeny.
    genus.list <- list()
    temp <- list()
    gu <- unique(genera[,1])
    for(i in 1:length(gu)){
        temp[[i]] <- genera[which(genera[,1] == gu[i]),]
        a <- sample(1:dim(temp[[i]])[1], 1)
        genus.list[[i]] <- temp[[i]][a,]
    }
    genus.list <- do.call(rbind,genus.list)
    
    ## Prune the tree only to keep those species:
    to.keep <- genus.list[,3]
    to.prune <- phy[[1]]$tip.label[!phy[[1]]$tip.label %in% to.keep]
    phy.genus <- lapply(1:100, FUN = function(x) drop.tip(phy[[x]], tip=to.prune))

    ## Change the tipname to the name of the genus:
    mm <- lapply(1:100, FUN = function(x) match(phy.genus[[x]]$tip.label, genus.list[,3]))
    for(i in 1:100){
        phy.genus[[i]]$tip.label <- genus.list[mm[[i]],1]
    }

    ## Setting all tree depths equal to 1.0
    phy.genus <- lapply(phy.genus, FUN = function(x) rescale(x, model = "depth", 1.0))

    return(phy.genus)
}
