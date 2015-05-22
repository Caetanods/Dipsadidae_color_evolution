my.medusa <-
function (x, time = TRUE, ...) 
{
    z <- x$zSummary
    shift <- list(cex = 1, bg = "gray", alpha = 0.75, col = "black", 
        lwd = 1)
    phy <- x$cache$phy
    mm <- match(phy$edge[, 2], z[, "dec"])
    edge.color <- z[mm, "partition"]
    plot(phy, edge.color = edge.color, ...)
    if (time) 
        axisPhylo()
    shifts <- z[(idx <- !is.na(z[, "shift"])), "dec"]
    if (length(shifts)) {
        ss <- z[idx, "shift"]
        ww <- character(nrow(phy$edge))
        ww[match(shifts, phy$edge[, 2])] <- ss
        xx <- numeric(nrow(phy$edge))
        xx[phy$edge[, 2] %in% shifts] <- 1
        edgelabels.auteur(NULL, frame = "circle", cex = ifelse(xx == 
            1, shift$cex, 1e-10), pch = ifelse(xx == 1, 21, NA), 
            bg = geiger:::.transparency(shift$bg, shift$alpha), col = shift$col, 
            lwd = shift$lwd)
        edgelabels.auteur(ww, frame = "none", cex = ifelse(xx == 
            1, shift$cex/3, 1e-10), pch = NA, col = shift$col)
    }
}
