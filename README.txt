README

The data and scripts presented here are an example of the analytical methods presented in "Changes in stream food web structure across a gradient of acid mine drainage increases local community stability" by Pomeranz et al. Ecology, in press. 

Archived scripts for all of the analysis and results presented in the paper are available from the main author upon request.



The full dataset presented in the article is available on datadryad, see main text for details. 


Data structure for this example
All of the data for this example is available in the example_data.RDS

The data is a list with 2 elements
1) community_data
data includes estimated dry weights (in grams) and local abundances (in number per m^2) for all taxa found from one site.
taxa = taxa name
dw.g = estimated dry weight in grams
density = number of individuals per meter squared
rel.ab = relative abundance of taxa
log10.dw = log base 10 of estimated dry weight

2) model_params - 
these are the estimated parameters for the trait matching model (See main text)

Script structure

1_example_script.R 
Loading necessary packages, functions and the example data.
Calculating the interaction probability matrices, including the Niche probabilities, pruning "forbidden" taxe, and accounting for neutral effects.
Estimate distribution of food web measures and stability metric by running n bernoulli trials. 

Example_functions.R
useful functions written for the analyses presented here

Stability_fxns_Sauve 
Supplementary R script with stability functions originally published in Sauve 2016
If using, please cite original publication:
Sauve, A. M. C., Thébault, E., Pocock, M. J. O., & Fontaine, C. (2016). How plants connect pollination and herbivory networks and their contribution to community stability. Ecology, 97(4), 908-917. doi.org/10.1890/15-0132.1
