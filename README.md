snake.colors.of.deception
=========================

R code for the article "Daniel S. Caetano, Laura R.V. Alencar, Paulo Passos, Marcio Martins, Felipe G. Grazziotin, & Hussam Zaher. The colors of deception: Evolution of warning signals in Neotropical snakes (Colubroidea: Dipsadidae)" To be submitted to Evolution.

This code will run the analysis, do simulations and make the graphs.

To run the complete analysis. Follow the order of scripts specified bellow:
 - 1) data-preparation.R (For all the alternative categorizations too. Then need to load the saved data in the BiSSE script.)
 - 2) BiSSE-analysis.R (This just run one of the BiSSE analysis. Need to collapse all estimation in this script. Thus the prepare data script need to prepare the dataset for all the runs and identify in a good way all the different types of analysis.)

The list of data used in all analyses:
 - "./data/beast_ingroup_100_posterior_sample.nex"
 - "./data/species_list.csv"
 - "./data/coloration_data.csv"
 - "./data/data_for_BiSSE.RData"

Note that only few changes are needed in order to run a similar analysis in your own dataset.
