## Does clades which show a higher proportion of contrasting species also show proportionaly higher richness?
## This is what you would expect from a model where traits are associated with speciation rates.
## Ok. But before I need to check if this expectation is true. Run a little bit of BiSSE simulations
##     and check if the thing make sense.

library(diversitree)
load("./data/data_for_BiSSE.RData")

state <- as.numeric(st[,2])
names(state) <- st[,1]

## State for the monotypic genera:
state <- state[!is.na(state)]
st.matrix <- data.frame(names(state), state, state, rep(1, length(state)), row.names = NULL)
colnames(st.matrix) <- colnames(unres)
st.matrix[which(st.matrix$n0 == 0), 2] <- 1
st.matrix[which(st.matrix$n1 == 1), 2] <- 0

## Create matrix with all the genera:
diver <- rbind(st.matrix, unres)

## Make regression:
prop <- diver$n1 / diver$Nc
plot(diver$Nc ~ prop, ylab = "Richness", xlab = "Proportion of contrasting spp.")
abline( lm(diver$Nc ~ prop ) )
summary(lm(diver$Nc ~ prop ))

## How this graph would look like with a single random draw from the same binomial
##     distribution I used to make the LRT simulations?

