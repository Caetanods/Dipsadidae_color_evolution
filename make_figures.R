## Script to make figures for the manuscript.
## Some of the figures were finished in InkScape for aesthetics. But the data points are
##     as ploted in R.

require(diversitree)

source("./functions/analysis.R")
load("./mcmc_BiSSE_results/results_100_phylo_bisse.RData")
load("./data/data_for_BiSSE.RData")

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
plot(density(a.p))

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

pdf(file = "Figure_3_parameter_estimates.pdf", width = 14, height = 7)
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
       , fill = c("#00000050","#ff000060","#00ff0040", "#0000ff40")
       , border = c("#00000070","#ff000080","#00ff0060", "#0000ff60")
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
