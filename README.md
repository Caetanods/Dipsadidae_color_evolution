snake.colors.of.deception
=========================

R code for the article "Daniel S. Caetano, Laura R.V. Alencar, Paulo Passos, Marcio Martins, Felipe G. Grazziotin, & Hussam Zaher. The colors of deception: Evolution of warning signals in Neotropical snakes (Colubroidea: Dipsadidae)" To be submitted to Evolution.

This code will run the analysis, do simulations and make the graphs.

To run the complete analysis. Follow the order of scripts specified bellow:
 - 1) data-preparation.R "Prepare all the data for the analyses. Including the different categorizations."
 - 2) BiSSE-analysis.R "Make the analysis for the main categorization and also the three other categories."

The list of data used in all analyses:
 - "./data/beast_ingroup_100_posterior_sample.nex"
 - "./data/species_list.csv"
 - "./data/coloration_data.csv"
 - "./data/data_for_BiSSE.RData"
 - "./data/data_for_BiSSE-alt.RData"

Note that only few changes are needed in order to run a similar analysis in your own dataset.
