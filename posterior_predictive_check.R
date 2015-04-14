## Perform the posterior predictive checks for the BiSSE analysis under the constrained and
##    unconstrained MCMC parameter estimates.

library(diversitree)

## Function that always return a tree with given parameters:
sims <- function(pars, taxa, x0){
    repeat
        {
            phy <- tree.bisse(pars, max.taxa = taxa, x0)
            if(!is.null(phy))
                {
                    break
                }
        }
    return(phy)
}

## Simulate data under the full model:
## The order of the parameters is lambda0; lambda1; mu0; mu1; q01; q10.
## The mean of the posteriors across all 100 sampled trees:
pars.full <- c(4.0, 7.4, 0.5, 0.7, 0.3, 3.0)
sim.full <- lapply(1:500, function(x) sims(pars = pars.full, taxa = 594, x0 = 1) )
sim.full.st1 <- sapply(1:500, function(x) sum(sim.full[[x]]$tip.state) )

## Note that with such parameters it is very unlikely to produce a tree
##    if the root state is set to 0. Uncomment the follow line to try.
## WARNING: It will run for a long time.
## sims(pars = pars.full, taxa = 594, x0 = 0)

## Now to simulate under the constrained model I just need to use the Mk estimates.
## Going to do the same as the previous BiSSE simulations. However, now the speciation and
##     extinction rates are going to be the same for both traits.
pars.cons <- c(5.5, 5.5, 0.5, 0.5, 0.2, 2.0)
sim.cons <- lapply(1:500, function(x) sims(pars = pars.cons, taxa = 594, x0 = 0) )
sims(pars = pars.cons, taxa = 594, x0 = 0)
sims(pars = pars.cons, taxa = 594, x0 = 1)
sim.cons.st1 <- sapply(1:500, function(x) sum(sim.cons[[x]]$tip.state) )

## Results.
## The empirical ratio is 213 / 594.
par(mfrow = c(1,2))
plot(density(( sim.cons.st1/594) *2 ))
plot(density(sim.full.st1/594))

par(mfrow = c(1,2))
hist(sim.cons.st1/594)
hist(sim.full.st1/594)

## It interesting that for both models the ratio is very different from the empirical.

