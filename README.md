snake.colors.of.deception
=========================

R code for the article "Daniel S. Caetano, Laura R.V. Alencar, Paulo Passos, Marcio Martins, Felipe G. Grazziotin, & Hussam Zaher. The colors of deception: Evolution of warning signals in Neotropical snakes (Colubroidea: Dipsadidae)" To be submitted to Evolution.

This code will run the analysis, do simulations and make the graphs.

To run the complete analysis. Follow the order of scripts specified bellow:
 - 1) data_preparation.R "Prepare all the data for the analyses. Including the different categorizations."
 - 2) BiSSE_analyses.R "Make the analysis for the main categorization and also the three other categories."
 - 3) model_selection_and_posterior_predictive_check.R "Make the model selection for all the MCMC results. Merge the posterior predictive analysis in this script too. May need to change the name for something more informative." "Need to make and include in the script the posterior predictive and DIC analysis for the alternative categorizations of the BiSSE analysis."
 - 4) medusa_analyses.R "Make the medusa analyses across the 100 phylogenies."
 - 5) BiSSE_null_model_simulations.R "Script with all the type I error simulations."
 - 6) make_figures.R "This need to reproduce all the figures in the last version of the manuscript."

The list of data used in all analyses:
 - "./data/beast_ingroup_100_posterior_sample.nex"
 - "./data/species_list.csv"
 - "./data/coloration_data.csv"
 - "./data/data_for_BiSSE.RData"
 - "./data/data_for_BiSSE-alt.RData"
 - "./data/dic_BiSSE_MCMC_results.RData"
 - "./data/asr_100_phylo.RData"
 - "./data/post_check.RData"

Note that only few changes are needed in order to run a similar analysis in your own dataset.
