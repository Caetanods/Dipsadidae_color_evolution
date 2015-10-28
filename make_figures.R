## Script to make figures for the manuscript.
## Some of the figures were finished in InkScape for aesthetics. But the data points are
##     as ploted in R.

library(diversitree)
library(geiger)

source("./functions/analysis.R")
load("./data/results_100_phylo_bisse.RData")
load("./data/data_for_BiSSE.RData")
tree.mcc <- read.tree("./data/ladder.tree.mcc.tre")

######################################################################################
## Figure 1 - Phylogeny, ancestral states and diversity at the tips.

## To make this figure need to first make ancestral state reconstruction on the MCC tree.

## Set depth of the MCC to 1 and solve polytomies for ancestral state estimation.
tree.mcc <- rescale(tree.mcc, model = "depth", 1.0)
tree.mcc <- multi2di(tree.mcc)

## Now get parameter estimates. For this going to use the grand mean of the mean estimates
##     from the posterior densities under the simple BiSSE model across the 100 phylogenies
##     sampled from the posterior distribution from Beast.
post <- seq(from=dim(mcmc.onerate[[1]][[1]])[1]/2, to=dim(mcmc.onerate[[1]][[1]])[1])
md <- list()
for(i in 1:length(mcmc.onerate)){
    md[[i]] <- sapply(mcmc.onerate[[i]][[1]][post,], mean)
}
md <- do.call(rbind, md)

## Simple plot to check the uncertainty in the parameter estimates across phylogenies:
boxplot(md[,-c(1,6)])

## Make the MLE estimate of the marginal ancestral reconstruction. This estimate use the simpler
##      BiSSE model that estimates a single net diversification for both traits. Thus, states do
##      not affect rates of extinction or speciation.

## In order to incorporate the variation in the estimate of the parameters for the simpler BiSSE
##      model in the ancestral state estimates we will make the MLE of the marginal ancestral state
##      under the 100 means from the different posterior distribution. The plot will show an
##      average of such estimates.

state <- as.numeric(st[,2])
names(state) <- st[,1]
lik <- make.bisse(tree.mcc, states = state, unresolved = unres)
lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
md <- md[,-c(1,6)]
st <- apply(md, 1, function(x) asr.marginal(lik, x)[2,])
st.avg2 <- rowMeans(st)
st.avg1 <- 1 - st.avg2
st <- unname(rbind(st.avg1, st.avg2))

## Making the figure:
## Preparing the bar with the proportion of states in each genera:
mm <- which(!tree.mcc$tip.label %in% unres[,1])
st.bar <- data.frame(tip.label = tree.mcc$tip.label[mm], state = state[tree.mcc$tip.label[mm]])
rownames(st.bar) <- NULL
st.bar <- data.frame(st.bar, n0 = rep(0, times = dim(st.bar)[1]), n1 = rep(0, times = dim(st.bar)[1]))
st.bar$n0[which(st.bar$state == 0)] <- 1
st.bar$n1[which(st.bar$state == 1)] <- 1
st.bar <- st.bar[,-2]

st.bar.total <- unres[,-4]
st.bar.total <- rbind(st.bar.total, st.bar)

mmm <- match(tree.mcc$tip.label, st.bar.total[,1])
tt <- st.bar.total[mmm,2] + st.bar.total[mmm,3]
names(tt) <- tree.mcc$tip.label

## Have aposematic and contrasting species.
st1 <- which(!st.bar.total[mmm,3] == 0 & !st.bar.total[mmm,2] == 0)
## Have only contrasting species.
st11 <- which(st.bar.total[mmm,2] == 0)
## Have only cryptic species.
st0 <- which(!st.bar.total[mmm,2] == 0)

pdf("Figure_1_Caetano_et_al_2016.pdf", width = 15, height = 20)
plot.phylo(tree.mcc, font = 4, label.offset = 0.03, cex = 1.0, adj = 1, x.lim = 2, edge.width = 3)
nodelabels(pie=t(st), piecol=c("grey","red"), cex=.35)
axisPhylo(side = 3, pos = c(70,0))

ad.t0 <- 1.275; incr <- 0.007; wd <- 0.2

## Add bars in three sections: All gray, red in top of gray, all red.
for(i in 1:length(st0)){ ## st0
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,2][st0[i]]), ad.t0+(incr*st.bar.total[mmm,2][st0[i]])
                , ad.t0), y = c(st0[i]-wd, st0[i]-wd, st0[i]+wd, st0[i]+wd), col = "grey", lend = 1)
}
for(i in 1:length(st1)){ ## st1 on top of st0
    polygon(x = c(ad.t0+(incr*st.bar.total[mmm,2][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]]))
          , y = c(st1[i]-wd, st1[i]-wd, st1[i]+wd, st1[i]+wd), col = "red", lend = 1)
}
for(i in 1:length(st11)){ ## st1 only
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,3][st11[i]]), ad.t0+(incr*st.bar.total[mmm,3][st11[i]])
                , ad.t0), y = c(st11[i]-wd, st11[i]-wd, st11[i]+wd, st11[i]+wd), col = "red", lend = 1)
}
dev.off()

######################################################################################
## Figure 2 - Not produced in R.

######################################################################################
## Figure 3 - Net diversification plot:

## Taking out a 50% burn-in and combining the posterior from the independent MCMC runs:
post.one.rate <- lapply(1:100, function(x) mcmc.onerate[[x]][[1]][5000:10000,] )
post.one.rate <- do.call(rbind, post.one.rate)
post.two.rate <- lapply(1:100, function(x) mcmc.tworate[[x]][[1]][5000:10000,] )
post.two.rate <- do.call(rbind, post.two.rate)

## Generate the prior distribution with same parameters used in the analyses:
start <- starting.point.bisse(tree.genus[[1]])
r <- 1 / 2 * (start[1] - start[3]) ## This is the 'rate' I used for the prior.
a <- seq(0,15, by = 0.01)
a.p <- r * exp(-r * a) ## Prior distribution for plots.

## Calculate the mean value of the distributions:
md.lambda0.one <- sapply(post.one.rate[2] - post.one.rate[3], median)
md.lambda0.two <- sapply(post.two.rate[2] - post.two.rate[4], median)
md.lambda1.two <- sapply(post.two.rate[3] - post.two.rate[5], median)
md.qs.one <- sapply(post.one.rate[4:5], median)
md.qs.two <- sapply(post.two.rate[6:7], median)

## Function to calculate the breaks:
to.break <- function(x, wd){
    y <- seq(floor(range(x))[1],range(x)[2], by = wd)
    bb <- length(y)
    if( y[bb] < max(x) ){
        step <- y[bb] - y[bb-1]
        y <- c(y, y[bb] + step)
    }
    return(y)
}

## Widths:
wd1 <- 0.15
wd2 <- 0.1

## Breaks for each parameter:
## First graph:
br.pr.1 <- to.break(a.p, wd1)
br.two.0 <- to.break(post.two.rate[,2] - post.two.rate[,4], wd1)
br.two.1 <- to.break(post.two.rate[,3] - post.two.rate[,5], wd1)
br.one.0 <- to.break(post.one.rate[,2] - post.one.rate[,3], wd1)
## Second graph:
br.pr.2 <- to.break(a.p, wd2)
br.two.q1 <- to.break(post.two.rate[,6], wd2)
br.two.q0 <- to.break(post.two.rate[,7], wd2)
br.one.q1 <- to.break(post.one.rate[,4], wd2)
br.one.q0 <- to.break(post.one.rate[,5], wd2)

pdf(file = "Figure_3_Caetano_et_al_2016.pdf", width = 14, height = 7)
par(mfrow = c(1,2))

## Make the graph with the net diversification rates (left):

h0 <- hist(a.p, breaks = br.pr.1, plot = FALSE)
h1 <- hist(post.one.rate[,2] - post.one.rate[,3], breaks = br.one.0, plot = FALSE)
h2 <- hist(post.two.rate[,2] - post.two.rate[,4], breaks = br.two.0, plot = FALSE)
h3 <- hist(post.two.rate[,3] - post.two.rate[,5], breaks = br.two.1, plot = FALSE)

h0$counts <- h0$counts/sum(h0$counts)
h1$counts <- h1$counts/sum(h1$counts)
h2$counts <- h2$counts/sum(h2$counts)
h3$counts <- h3$counts/sum(h3$counts)

## Blue (prior):
plot(h0, main = "", xlab = "", ylab = ""
     , col = "#0000ff40", border = "#0000ff60", xlim = c(0,12), ylim = c(0,1.0), axes = FALSE)
## Yellow:
plot(h1, xlab = "", ylab = "", main = "", add = TRUE, col = "#fecd0060", border = "#fecd0080")
## Gray:
plot(h2, xlab = "", ylab = "", main = "", add = TRUE, col = "#00000040", border = "#00000060")
## Red:
plot(h3, xlab = "", ylab = "", main = "", add = TRUE, col = "#ff000060", border = "#ff000080")

## The 95% CI lines:
div1 <- quantile(ecdf(post.two.rate[,2] - post.two.rate[,4]), probs = c(0.025, 0.975))
div2 <- quantile(ecdf(post.two.rate[,3] - post.two.rate[,5]), probs = c(0.025, 0.975))
div3 <- quantile(ecdf(post.one.rate[,2] - post.one.rate[,3]), probs = c(0.025, 0.975))

segments(x0 = div1[1], y0 = -0.015, x1 = div1[2], y1 = -0.015, col = "#00000060", lwd = 4, lend = 1)
segments(x0 = div2[1], y0 = -0.015, x1 = div2[2], y1 = -0.015, col = "#ff000080", lwd = 4, lend = 1)
segments(x0 = div3[1], y0 = -0.015-0.018, x1 = div3[2], y1 = -0.015-0.018, col = "#fecd0060", lwd = 4, lend = 1)

## The mean lines:
segments(x0 = md.lambda0.two, y0 = 0, x1 = md.lambda0.two, y1 = 0.25, col = "#00000060"
         , lwd = 2.5, lend = 1, lty = 2)
segments(x0 = md.lambda1.two, y0 = 0, x1 = md.lambda1.two, y1 = 0.25, col = "#ff000080"
         , lwd = 2.5, lend = 1, lty = 2)
segments(x0 = md.lambda0.one, y0 = 0, x1 = md.lambda0.one, y1 = 0.25, col = "#fecd0060"
         , lwd = 2.5, lend = 1, lty = 2)
## Median values in the top of the lines:
tx <- c(md.lambda0.two, md.lambda1.two, md.lambda0.one)
text(x = tx, y = 0.27, labels = as.character(round(tx, digits = 1)))

axis(side = 2); axis(side = 1, pos = -0.05)
mtext("Net diversification", 1, line=2.5, cex = 1.5)
mtext("Probability density", 2, line=2.5, cex = 1.5)
legend(x = 6.0, y = 0.8,legend = c("cryptic","contrasting","color independent", "prior")
       , fill = c("#00000050","#ff000060","#fecd0060", "#0000ff40")
       , border = c("#00000070","#ff000080","#fecd0080", "#0000ff60")
       , bty = "n", ncol = 1, cex = 1.5)

## Make the graph with the transition rates (right):

h1 <- hist(a.p, breaks = br.pr.2, plot = FALSE)
h2 <- hist(post.two.rate[,6], breaks = br.two.q1, plot = FALSE)
h3 <- hist(post.two.rate[,7], breaks = br.two.q0, plot = FALSE)
h4 <- hist(post.one.rate[,4], breaks = br.one.q1, plot = FALSE)
h5 <- hist(post.one.rate[,5], breaks = br.one.q0, plot = FALSE)

h1$counts <- h1$counts/sum(h1$counts)
h2$counts <- h2$counts/sum(h2$counts)
h3$counts <- h3$counts/sum(h3$counts)
h4$counts <- h4$counts/sum(h4$counts)
h5$counts <- h5$counts/sum(h5$counts)

## Blue (prior):
plot(h1, main = "", xlab = "", ylab = "", col = "#0000ff40", border = "#0000ff60"
   , xlim = c(0,8), ylim = c(0,1), axes = FALSE)
## Light red:
plot(h2, xlab = "", ylab = "", add = TRUE, col = "#ff000040", border = "#ff000060", axes = FALSE)
## Light gray:
plot(h3 , xlab = "", ylab = "", add = TRUE, col = "#00000040", border = "#00000060", axes = FALSE)
## Dark red:
plot(h4, xlab = "", ylab = "", add = TRUE, col = "#ff000070", border = "#ff000090", axes = FALSE)
## Dark grey:
plot(h5 , xlab = "", ylab = "", add = TRUE, col = "#00000070", border = "#00000090", axes = FALSE)

tra1 <- quantile(ecdf(post.two.rate[,6]), probs = c(0.025, 0.975))
tra2 <- quantile(ecdf(post.two.rate[,7]), probs = c(0.025, 0.975))
tra3 <- quantile(ecdf(post.one.rate[,4]), probs = c(0.025, 0.975))
tra4 <- quantile(ecdf(post.one.rate[,5]), probs = c(0.025, 0.975))

segments(x0 = tra1[1], y0 = -0.015, x1 = tra1[2], y1 = -0.015, col = "#ff000060", lwd = 4, lend = 1)
segments(x0 = tra2[1], y0 = -0.015, x1 = tra2[2], y1 = -0.015, col = "#00000060", lwd = 4, lend = 1)
segments(x0 = tra3[1], y0 = -0.015-0.018, x1 = tra3[2], y1 = -0.015-0.018, col = "#ff000090", lwd = 4, lend = 1)
segments(x0 = tra4[1], y0 = -0.015-0.018, x1 = tra4[2], y1 = -0.015-0.018, col = "#00000090", lwd = 4, lend = 1)

segments(x0 = md.qs.two, y0 = 0, x1 = md.qs.two, y1 = c(0.4, 0.25), col = c("#ff000060","#00000060")
         , lwd = 2.5, lend = 1, lty = 2)
segments(x0 = md.qs.one, y0 = 0, x1 = md.qs.one, y1 = c(0.4, 0.25), col = c("#ff000090","#00000090")
         , lwd = 2.5, lend = 1, lty = 2)
## Little adjustment to the text of "cryptic to aposematic".
text(x = c(md.qs.two[1]+0.1, md.qs.two[2]), y = c(0.42, 0.27), labels = as.character(round(md.qs.two, digits = 1)))
text(x = c(md.qs.one[1]+0.1, md.qs.one[2]), y = c(0.42, 0.27), labels = as.character(round(md.qs.one, digits = 1)))

axis(side = 2); axis(side = 1, pos = -0.05)
mtext("State transition rates", 1, line=2.5, cex = 1.5)
mtext("Probability density", 2, line=2.5, cex = 1.5)
legend(x = 3.8, y = 0.8,legend = c("cryptic to contrasting (q_01)","contrasting to cryptic (q_10)", "prior")
       , fill = c("#ff000060","#00000050","#0000ff40")
       , border = c("#ff000080","#00000070","#0000ff60")
       , bty = "n", ncol = 1, cex = 1.5)

dev.off()
######################################################################################
## FIGURES FOR THE SUPPLEMENTARY INFORMATION

## Phylogeny, ancestral states and diversity at the tips using the mimics vs. non-mimics category.
## Note that this figure is just a version of the Figure 1 included in the main test.

## Get the data and parameters for the altC analysis.
load("./data/results_bisse_mcmc_alt_ABC.RData")
load("./data/data_for_BiSSE-alt.RData")

## To make this figure need to first make ancestral state reconstruction on the MCC tree.
## Set depth of the MCC to 1 and solve polytomies for ancestral state estimation.
tree.mcc <- rescale(tree.mcc, model = "depth", 1.0)
tree.mcc <- multi2di(tree.mcc)

## Now get parameter estimates. For this going to use the grand mean of the mean estimates
##     from the posterior densities under the simple BiSSE model across the 100 phylogenies
##     sampled from the posterior distribution from Beast.
post <- seq(from=dim(mcmc.onerate.C[[1]])[1]/2, to=dim(mcmc.onerate.C[[1]])[1])
md <- list()
for(i in 1:length(mcmc.onerate.C)){
    md[[i]] <- sapply(mcmc.onerate.C[[i]][post,], mean)
}
md <- do.call(rbind, md)

## Simple plot to check the uncertainty in the parameter estimates across phylogenies:
boxplot(md[,-c(1,2,7)])

md

stateC <- as.numeric(st.alt[[3]][,2])
names(stateC) <- st.alt[[3]][,1]

lik <- make.bisse(tree.mcc, states = stateC, unresolved = unres.alt[[3]])
lik <- constrain(lik, lambda1~lambda0, mu1~mu0)
md <- md[,-c(1,2,7)]
st <- apply(md, 1, function(x) asr.marginal(lik, x)[2,])
st.avg2 <- rowMeans(st)
st.avg1 <- 1 - st.avg2
st <- unname(rbind(st.avg1, st.avg2))

## Making the figure:
## Preparing the bar with the proportion of states in each genera:
mm <- which(!tree.mcc$tip.label %in% unres.alt[[3]][,1])
st.bar <- data.frame(tip.label = tree.mcc$tip.label[mm], state = stateC[tree.mcc$tip.label[mm]])
rownames(st.bar) <- NULL
st.bar <- data.frame(st.bar, n0 = rep(0, times = dim(st.bar)[1]), n1 = rep(0, times = dim(st.bar)[1]))
st.bar$n0[which(st.bar$state == 0)] <- 1
st.bar$n1[which(st.bar$state == 1)] <- 1
st.bar <- st.bar[,-2]

st.bar.total <- unres.alt[[3]][,-4]
st.bar.total <- rbind(st.bar.total, st.bar)

mmm <- match(tree.mcc$tip.label, st.bar.total[,1])
tt <- st.bar.total[mmm,2] + st.bar.total[mmm,3]
names(tt) <- tree.mcc$tip.label

## Have mimic and non-mimic species.
st1 <- which(!st.bar.total[mmm,3] == 0 & !st.bar.total[mmm,2] == 0)
## Have only mimic species.
st11 <- which(st.bar.total[mmm,2] == 0)
## Have only non-mimic species.
st0 <- which(!st.bar.total[mmm,2] == 0)

pdf("Figure_SI_Ancestral_state_mimic_non-mimic_Caetano_et_al_2016.pdf", width = 15, height = 20)
plot.phylo(tree.mcc, font = 4, label.offset = 0.03, cex = 1.0, adj = 1, x.lim = 2, edge.width = 3)
nodelabels(pie=t(st), piecol=c("grey","red"), cex=.35)
axisPhylo(side = 3, pos = c(70,0))

ad.t0 <- 1.275; incr <- 0.007; wd <- 0.2

## Add bars in three sections: All gray, red in top of gray, all red.
for(i in 1:length(st0)){ ## st0
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,2][st0[i]]), ad.t0+(incr*st.bar.total[mmm,2][st0[i]])
                , ad.t0), y = c(st0[i]-wd, st0[i]-wd, st0[i]+wd, st0[i]+wd), col = "grey", lend = 1)
}
for(i in 1:length(st1)){ ## st1 on top of st0
    polygon(x = c(ad.t0+(incr*st.bar.total[mmm,2][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]])+(incr*st.bar.total[mmm,3][st1[i]])
                , ad.t0+(incr*st.bar.total[mmm,2][st1[i]]))
          , y = c(st1[i]-wd, st1[i]-wd, st1[i]+wd, st1[i]+wd), col = "red", lend = 1)
}
for(i in 1:length(st11)){ ## st1 only
    polygon(x = c(ad.t0, ad.t0+(incr*st.bar.total[mmm,3][st11[i]]), ad.t0+(incr*st.bar.total[mmm,3][st11[i]])
                , ad.t0), y = c(st11[i]-wd, st11[i]-wd, st11[i]+wd, st11[i]+wd), col = "red", lend = 1)
}
dev.off()
