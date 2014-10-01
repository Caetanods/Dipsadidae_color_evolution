## This script will run the alternative categorization for the colors of the Dipsadidae snakes.
## There are two alternative categorizations:
## a) Species which show ventral contrasting colors as "CR".
## b) Species which only the juvenile show contrastiong colors as "CR".

library(diversitree)

## Load the data and the functions:
source("./functions/data-prepare.R")
source("./functions/analysis.R")
data <- read.csv("./data/coloration_data.csv", as.is = TRUE)[,-c(1,4)]
load("./data/data_for_BiSSE.RData")
rm(st, unres)

###########################
## Alternative categorization A: VEN -> CR

altA <- data$Category

## Tranforming the polymorphisms into single categories.
## Change those transformations for alternative codding.

## Aposematic category:
altA[which(altA == "CON")] <- 1

## Changes: VEN -> CR:
altA[which(altA == "VEN")] <- 0
altA[which(altA == "CON+VEN")] <- 0

## Cryptic category:
altA[which(altA == "CR")] <- 0
altA[which(altA == "VIP")] <- 0
altA[which(altA == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorfic here.
altA[which(!altA == "1")] <- 0
data$altA <- altA
###########################

###########################
## Alternative categorization B: juvenile only -> CR

altB <- data$Category

## Tranforming the polymorphisms into single categories.
## Change those transformations for alternative codding.

## Aposematic category:
altA[which(altA == "CON")] <- 1
altA[which(altA == "VEN")] <- 1
altA[which(altA == "CON+VEN")] <- 1

## Changes: VEN -> CR:

## Cryptic category:
altA[which(altA == "CR")] <- 0
altA[which(altA == "VIP")] <- 0
altA[which(altA == "CR+VIP")] <- 0

## Since there are cryptic forms among the polymorfic here.
altA[which(!altA == "1")] <- 0
data$altA <- altA
