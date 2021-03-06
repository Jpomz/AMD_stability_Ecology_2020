# 1_example_script

# this shows an example of how to predict interaction probabilities based on body size and relative abundance. 
# This analysis utilizes adjacency matrix (A) structure to represent interactions. Matrices are ordered by increasing average body size from left to right and top to bottom. Columns represent species in their role as consumers, and rows represent their role as resources. 
# i.e., Aij = 1 when column j consumes row i, and Aij = 0 if column j does not consumer row i. See methods in main text.

# Load functions ----------------------------------------------------------

# load predict.niche.prob() function from Bartomeus et al. 2016
source("functions/predict.niche.prob.R")

# stability functions from Sauve et al. 2016
source("functions/stability_fxns_Sauve.R")

# load functions written for this manuscript
source("functions/Example_functions.R")

# Load data ---------------------------------------------------------------
example_data <- readRDS("data/example_data.RDS")

# make sure data is sorted by increasing body size
example_data$community_data


# Interaction probability - NICHE -----------------------------------------

# make matrix with all possible pairwise combinations of log10 taxa dryweights
M <- expand.grid(example_data$community_data$log10.dw,
                 example_data$community_data$log10.dw)

# estimate probability of interaction based on body sizes
# This function solves equation 5 from Bartomeus et al. 2016, and predicts the probabiility of an interaction between two new species based on their body sizes
# if using, please cite Bartomeus et al. 2016
# the "model_params" are from a model parameterized using the traitmatch package. 
# The model was parameterized using individual feeding data from Broadstone Stream and Tadnoll Brook (Woodward et al. 2010). Data was provided by Guy Woodward and Iwan Jones (personal communication, contact these authors for access to data)
link.probs <- predict.niche.prob(pars = example_data$model_params,
                                 M[[1]],
                                 M[[2]],
                                 replicates = 1)[[1]]

# this is a vector of all of the probabilities of pairwise interactions. We are going to convert this to a matrix. 
prob.matr <- matrix(link.probs,
                    nrow = sqrt(length(link.probs)),
                    ncol = sqrt(length(link.probs)))

# add taxa names to matrix rows and columns 
# this is necessary for pruning niche forbidden links below
dimnames(prob.matr) <- list(example_data$community_data$taxa,
                            example_data$community_data$taxa)
# print out sample of matrix
prob.matr[1:5, 1:5]


# Relative abundance - NEUTRAL EFFECTS ------------------------------------

# Now we are going to make a matrix of relative abundances N, where Nij = Relative abundance of species i times relative abundance of species j 
# See the "get_rel_ab()" function in "MS_functions.R" for more details
rel.ab.matr <- get_rel_ab(vec = example_data$community_data$rel.ab,
                          taxa = example_data$community_data$taxa)

# sample of relative abundance matrix
rel.ab.matr[1:5, 1:5]

# rescale the relative abundances in this matrix to be from 0.5 to 1
# this means that the most abundant species pair are likely to interact, while the rarest pairs are less likely to interact. 
# See the "scalexy()" function in "MS_functions.R" for more details
rel.ab.matr <- scalexy(rel.ab.matr, min = 0.5, max = 1)
# sample of rescaled relative abundance matrix
rel.ab.matr[1:5, 1:5]


# Final probability matrices ----------------------------------------------

# 1) prune niche forbidden links e.g. Pomeranz et al. 2019 Methods Eco Evo
# 2) estimate link probability by multiplying NICHE probability (link.probs) and NEUTRAL probability (rel.ab.matr)

# 1) prune niche "forbidden" taxa based on morphology
# e.g. "brushing/scraping" mouthparts 
taxa.forbid <- c("Acari", "Austrosimulium", "Coloburiscus",
                 "Deleatidium", "Elmidae", "Helicopsyche",
                 "Hydraenidae","Nesameletus", "Oligochaeta",
                 "Olinga", "Oxyethira", "Potamopyrgus",
                 "Spaniocerca", "Zelandobius")

# set column values to 0 for forbidden taxa (i.e. these taxa are non-predatory, so set their consumer roles to 0)
# see rm_niche() function in MS_functions.R
# make sure your probability matrix has taxa in column names, otherwise it won't do wnything. 
prob.matr.pruned <- rm_niche(prob.matr,
                             f.taxa = taxa.forbid)

# Pruned matrix, note that some columns are all 0's now
prob.matr.pruned[1:5, 1:5]

# 2) NEUTRAL EFFECTS
# "weighting" interaction probabilities based on relative abundances
# e.g. more abundant pairs are more likely to interact
prob.matr.neutral <- prob.matr.pruned * rel.ab.matr

# print out sample of 
prob.matr.neutral[1:5, 1:5]

# rescale probabilities to be between 0.01 and 0.99
prob.matr.final <- scalexy(prob.matr.neutral,
                             min = 0.01, 
                             max = 0.99)

# final probability matrix
prob.matr.final[1:5, 1:5]

# Adjacency matrices ------------------------------------------------------

# The b_trial() function in Example_functions.R takes a square, probability matrix, and returns a binary, adjacency matrix based on those probabilities

b_trial(prob.matr.final)

# you can save these iterations as objects, and then get standard food web measures from them. 
A1 <- b_trial(prob.matr.final)
Get.web.stats(A1)

# you can also estimate stability from them, shown below. 

# Jacobian matrices and stability -----------------------------------------

# the "get_measures()" function in the MS_functions.R script takes a probability matrix, and runs numerous "trials". For each trial, it returns the estimated stability metric "stab", as well as a number of standard food web measures. 

# this analsis and functions relies heavily on the stability functions of Sauve et al. 2016. If using, please cite the original publication: Sauve, A. M. C., Th�bault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1

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

# references

# Bartomeus I., Gravel D., Tylianakis J.M., Aizen M.A., Dickie I.A. & Bernard-Verdier M. (2016). A common framework for identifying linkage rules across different types of interactions. Functional Ecology 30, 1894-1903. https://doi.org/10.1111/1365-2435.12666

# Pomeranz J.P.F., Thompson R.M., Poisot T. & Harding J.S. (2019). Inferring predator-prey interactions in food webs. Methods in Ecology and Evolution 10, 356-367. https://doi.org/10.1111/2041-210X.13125

# Sauve A.M.C., Th�bault E., Pocock M.J.O. & Fontaine C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology 97, 908-917 doi.org/10.1890/15-0132.1

# Woodward G., Blanchard J.L., Lauridsen R.B., Edwards F.K., Jones J.I., Figueroa D.H., et al. (2010). Individual-Based Food Webs.: Species Identity , Body Size and Sampling Effects Individual-Based Food Webs.: Species Identity , Body Size and Sampling Effects. Advances In Ecological Research 43, 211-266. https://doi.org/10.1016/B978-0-12-385005-8.00006-X
