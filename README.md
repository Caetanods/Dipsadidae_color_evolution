snake.colors.of.deception
=========================

R code for the article "Daniel S. Caetano, Laura R.V. Alencar, Paulo Passos, Marcio Martins, Felipe G. Grazziotin, & Hussam Zaher. The colors of deception: Evolution of warning signals in Neotropical snakes (Colubroidea: Dipsadidae)" To be submitted to Evolution.

Here are all the data and scripts necessary to reproduce all the analyses and figures. Each script is self explanatory. However the following order is recommended for a better understanding of the steps applied to the analyses.

Note that some of the analyses and simulations take a long time to run. In those cases we present the code in a commented form immediately preceeded by a warning. The user can easily skip blocks of commented code with such analyses and load the results instead of running the simulations again.

Run the scripts in the following order:

 - 1) data_preparation.R
 Prepare all the data for the analyses.
 - 2) BiSSE_analyses.R
 Make all analyses. Including the main categorization and the three alternative definitions for the coloration patterns.
 - 3) model_selection_and_posterior_predictive_check.R
 Make the model selection for all the MCMC results using Deviance Information Criteria (DIC). Also make posterior predictive checks to test the adequacy of the model to the empirical data.
 - 4) medusa_analyses.R
 Perform MEDUSA analyses across the 100 phylogenies sampled from the posterior distribution of a Beast 1.8 phylogeny estimation with a concatenated dataset.
 - 5) BiSSE_null_model_simulations.R
 Perform Likelihood Ratio Test (LRT) simulations under the BiSSE null model following Rabosky and Goldberg, 2015.
 - 6) make_figures.R
 Produce all the figures reported in the manucript.
