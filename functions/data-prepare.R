to.genus.tree <- function(phy, genera){
    ## Picking one species per genus randomly:
    ## phy = multiphylo object
    ## genera = matrix or data.frame. first column with the genera names (or any unique identifier).
    ##        Second column with the tips names, this need to be the same name of the tips of the
    ##        original phylogeny.
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
    to.keep <- genus.list[,2]
    to.prune <- phy[[1]]$tip.label[!phy[[1]]$tip.label %in% to.keep]
    phy.genus <- lapply(1:100, FUN = function(x) drop.tip(phy[[x]], tip=to.prune))

    ## Change the tipname to the name of the genus:
    mm <- lapply(1:100, FUN = function(x) match(phy.genus[[x]]$tip.label, genus.list[,2]))
    for(i in 1:100){
        phy.genus[[i]]$tip.label <- as.character(genus.list[mm[[i]],1])
    }

    ## Setting all tree depths equal to 1.0
    phy.genus <- lapply(phy.genus, FUN = function(x) rescale(x, model = "depth", 1.0))

    return(phy.genus)
}

tree.pruner <- function(phy, outgroup){
    label <- phy$tip.label
    outgroup[!outgroup %in% label]
    prune <- which(label %in% outgroup)
    ingroup <- label[-prune]
    if(!is.monophyletic(phy, tips = ingroup)) stop("ingroup is not monophyletic")
    cphy <- drop.tip(phy, tip = outgroup)
    return(cphy)
}

## Below are functions special to the "no_effect_model_test.R" script.
make.matrix <- function(uu, st){
    mm <- data.frame(uu$Var1, uu$Freq - st, st, uu$Freq, stringsAsFactors = FALSE)
    colnames(mm) <- c("tip.label","n0","n1","Nc")
    rownames(mm) <- NULL
    return(mm)
}

make.unres <- function(null.states){
    genus.table <- data.frame(table(null.states[,1]))
    i <- sapply(genus.table, is.factor)
    genus.table[i] <- lapply(genus.table[i], as.character)

    one.species <- genus.table[which(genus.table[,2] == 1),1]
    species.table <- null.states[which(!null.states[,1] %in% one.species),]
    
    uu <- data.frame(table(as.character(species.table[,1])))
    i <- sapply(uu, is.factor)
    uu[i] <- lapply(uu[i], as.character)

    st.1 <- t(sapply(uu$Var1, function(x)
        colSums(species.table[which(species.table[,1] == x),-c(1,2)]) ))

    unres <- lapply(1:ncol(st.1), function(x) make.matrix(uu, st.1[,x]) )
    return(unres)
}

make.bisse.states <- function(null.states){
    genus.table <- data.frame(table(null.states[,1]))
    i <- sapply(genus.table, is.factor)
    genus.table[i] <- lapply(genus.table[i], as.character)

    one.species <- genus.table[which(genus.table[,2] == 1),1]
    other <- genus.table[which(!genus.table[,2] == 1),1]
    one.states <- null.states[which(null.states[,1] %in% one.species),]

    st <- lapply(3:ncol(one.states), function(x) one.states[,c(1,x)] )
    for(i in 1:length(st)) colnames(st[[i]]) <- c("Genus","BISSE")
    app <- cbind(other, NA)
    colnames(app) <- c("Genus","BISSE")
    st.bisse <- lapply(1:length(st), function(x) as.matrix(rbind(st[[x]], app)) )
    states <- list()
    for(i in 1:length(st)){
        states[[i]] <- as.numeric(st.bisse[[i]][,2])
        names(states[[i]]) <- unname(st.bisse[[i]][,1])
    }

    return(states)
}
