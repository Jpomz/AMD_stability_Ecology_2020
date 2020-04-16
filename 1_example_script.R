# 1_example_script


# Load libraries ----------------------------------------------------------

library(dplyr)
library(GenSA)
library(remotes)
# #SDMTools no longer on CRAN - must install archived version as below
# install.packages("https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz", repo=NULL, type="source")
# library(SDMTools)
# 
# #install "traitmatch" package - https://github.com/ibartomeus/traitmatch - WILL NOT INSTALL UNLESS SDMTools is installed first
# remotes::install_github("ibartomeus/traitmatch", dependencies = T)
library(traitmatch)

# load predict.niche.prob() function
source("functions/predict.niche.prob.R")

# stability functions from Sauve et al. 2016
source("functions/stability_fxns_Sauve.R")

# load functions written for this manuscript
source("functions/Example_functions.R")

# Load data ---------------------------------------------------------------
example_data <- readRDS("data/example_data.RDS")


# Interaction probability - NICHE -----------------------------------------

# make matrix with all possible pairwise combinations of taxa dryweights
M <- expand.grid(example_data$community_data$log10.dw,
                 example_data$community_data$log10.dw)

# estimate probability of interaction based on body sizes
# this uses the traitmatch package, see Bartomeus et al. 2016 for more details. 
link.probs <- predict.niche.prob(pars = example_data$model_params,
                                 M[[1]],
                                 M[[2]],
                                 replicates = 1)[[1]]

# this is a vector of all of the probabilities of pairwise interactions. We are going to convert this to a matrix. 
prob.matr <- matrix(link.probs,
                    nrow = sqrt(length(link.probs)),
                    ncol = sqrt(length(link.probs)))



# Relative abundance - NEUTRAL EFFECTS ------------------------------------

# See the "get_rel_ab()" function in "MS_functions.R" for more details
rel.ab.matr <- get_rel_ab(vec = example_data$community_data$rel.ab,
                          taxa = example_data$community_data$taxa)


# rescale the relative abundances in this matrix to be from 0.5 to 1
# See the "scalexy()" function in "MS_functions.R" for more details
rel.ab.matr <- scalexy(rel.ab.matr, min = 0.5, max = 1)


# Final probability matrices ----------------------------------------------

# 1) prune niche forbidden links e.g. Pomeranz et al. 2019 Methods Eco Evo
# 2) estimate link probability by multiplying NICHE probability (link.probs) and NEUTRAL probability (rel.ab.matr)

# 1) prune niche "forbidden" taxa based on morphology
# e.g. "brushing/scraping" mouthparts, 
taxa.forbid <- c("Acari", "Austrosimulium", "Coloburiscus",
                 "Deleatidium", "Elmidae", "Helicopsyche",
                 "Hydraenidae","Nesameletus", "Oligochaeta",
                 "Olinga", "Oxyethira", "Potamopyrgus",
                 "Spaniocerca", "Zelandobius")

# set column values to 0 for forbidden taxa
# see rm_niche() function in MS_functions.R
prob.matr.pruned <- rm_niche(prob.matr,
                             taxa = taxa.forbid)

# 2) NEUTRAL EFFECTS
# "weighting" interaction probabilities based on relative abundances
# e.g. more abundant pairs are more likely to interact
prob.matr.neutral <- prob.matr.pruned * rel.ab.matr

# rescale probabilities to be between 0.01 and 0.99
prob.matr.final <- scalexy(prob.matr.neutral,
                             min = 0.01, 
                             max = 0.99)


# Adjacency matrices ------------------------------------------------------

# The b_trial() function in Example_functions.R takes a square, probability matrix, and returns a binary, adjacency matrix based on those probabilities

b_trial(prob.matr.final)

# you can save these iterations as objects, and then get standard food web measures from them. 
A1 <- b_trial(prob.matr.final)
Get.web.stats(A1)

# you can also estimate stability from them, shown below. 

# Jacobian matrices and stability -----------------------------------------

# the "get_measures()" function in the MS_functions.R script takes a probability matrix, and runs numerous "trials". For each trial, it returns the estimated stability metric "stab", as well as a number of standard food web measures. 

# this analsis and functions relies heavily on the stability functions of Sauve et al. 2016. If using, please cite the original publication: Sauve, A. M. C., Thébault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1

# The example below is based on random interaction strengths 

# the function allows you to "scale" interaction strengths based on relative body size
# i.e. argument "scaled.Jij = TRUE"
# It can also correlate positive and negative interaction strengths  
# i.e. "correlate.Jij = TRUE"
# default "correlate.value" = 0.7, i.e. positive interaction = 0.7 * negative interaction, but you can change this setting. 
# You can also set both scaled and correlate to TRUE to examine effects of both

# You can control the number of replications with the "trials" argument. 

random.J <- get_measures(prob.matr.final,
                         s2 = 2,
                         trials = 10,
                         scale.Jij = FALSE,
                         correlate.Jij = TRUE)

random.J
# the output is a list of length = trials
# each element in the list contains the following information:
# stab = stability metric from Sauve et al. 2016 (lower values are "more" stable)
# S = number of species (taxa) in web, these should be identical for each run per site
# L = number of links in that iteration - this will change for each iteration due to the random nature of bernoulli trials
# C = connectance = L / S^2
# B = number of basal taxa (i.e., prey, not predators)
# I = number of intermediate taxa (i.e., both prey and predators)
# T = number of top taxa (i.e., predators, but not prey)
