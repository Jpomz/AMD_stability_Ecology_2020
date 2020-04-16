# 1_example_script

# Load packages for analysis
# you may need to install the packages first
# to do this, run "install.packages("package_name"), then load it
# for example, first run: install.packages("dplyr")
# then run: library(dplyr)
#library(dplyr)

# to install "traitmatch" package, see:
# https://github.com/ibartomeus/traitmatch
# note that the syntax of the command:
# install_github("traitmatch", "ibartomeus")
# is deprecated. 
# replace with :
# install_github("ibartomeus/traitmatch")
library(traitmatch)

# load predict.niche.prob() function
source("functions/predict.niche.prob.R")

# stability functions from Sauve et al. 2016
source("functions/stability_fxns_Sauve.R")

# load functions written for this manuscript
source("functions/Example_functions.R")

# load data
#saveRDS(example_data, "data/example_data.RDS")
example_data <- readRDS("data/example_data.RDS")

# probability of interactions body size ####

# make matrix with all possible pairwise combinations of taxa dryweights
M <- expand.grid(example_data$community_data$log10.dw, example_data$community_data$log10.dw)

# estimate probability of interaction based on body sizes
# this uses the traitmatch package, see Bartomeus et al. 2016 for more details. 
link.probs <- predict.niche.prob(pars = example_data$model_params,
                                 M[[1]],
                                 M[[2]],
                                 replicates = 1)[[1]]

# this is a vector of all of the probabilities of pairwise interactions. We are going to convert this to a matrix. 
prob.matr <- matrix(link.probs, sqrt(length(link.probs)), sqrt(length(link.probs)))


# relative abundance matrix ####
# See the "get_rel_ab()" function in "MS_functions.R" for more details
rel.ab.matr <- get_rel_ab(vec = example_data$community_data$rel.ab,
                          taxa = example_data$community_data$taxa)


# rescale the relative abundances in this matrix to be from 0.5 to 1
# See the "scalexy()" function in "MS_functions.R" for more details
rel.ab.matr <- scalexy(rel.ab.matr, min = 0.5, max = 1)

# final probability matrices ####
# 1) prune niche forbidden links e.g. Pomeranz et al. 2019 Methods Eco Evo
# 2) estimate link probability by multiplying niche probability (link.probs) and neutral probability (rel.ab.matr)

# 1) prune niche "forbidden" taxa based on morphology
# e.g. "brushing/scraping" mouthparts, 
taxa.forbid <- c("Acari", "Austrosimulium", "Coloburiscus", "Deleatidium", "Elmidae", "Helicopsyche", "Hydraenidae", "Nesameletus","Oligochaeta", "Olinga", "Oxyethira", "Potamopyrgus", "Spaniocerca", "Zelandobius")

# set column values to 0 for forbidden taxa
# see rm_niche() function in MS_functions.R
prob.matr.pruned <- rm_niche(prob.matr, taxa = taxa.forbid)

# 2) neutral effects
# "weighting" interaction probabilities based on relative abundances
# e.g. more abundant pairs are more likely to interact
prob.matr.neutral <- prob.matr.pruned * rel.ab.matr

# rescale probabilities to be between 0.01 and 0.99
prob.matr.final <- scalexy(prob.matr.neutral,
                             min = 0.01, 
                             max = 0.99)

# Jacobian matrices and stability ####
# these use the "get_measures()" function in the MS_functions.R script. 

# this allows relies heavily on the stability functions of Sauve et al. 2016. If using, please cite the original publication: Sauve, A. M. C., Thébault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1

# This example shows random interaction strengths 

# the function allows you to "scale" interaction strengths based on relative body size
# i.e. argument "scaled.Jij = TRUE"
# It can also correlate positive and negative interaction strengths  
# i.e. "correlate.Jij = TRUE"
# default "correlate.value" = 0.7, i.e. positive interaction = 0.7 * negative interaction
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
