# Example functions 
# jfpomeranz@gmail.com
# June 2018

# These are useful functions for running the analyses presented in Pomeranz et al. 2020 Ecology

# bernoulli trials
# this function takes a square, probability matrix, and returns a binary adjacency matrix of the same dimensions
b_trial <- function (prob_matr){
  sum.out = 0
  # make sure that all matrices have at least one link
  while(sum.out == 0){
    out <- matrix(rbinom(length(prob_matr),
                         1,
                         prob = prob_matr),
                  ncol = ncol(prob_matr),
                  nrow = nrow(prob_matr))
    sum.out = sum(out[upper.tri(out)])
  }
  out
}

# function to calculate relative abundance matrices
# this function takes a vector of abundance data and a character vector of taxa names. It returns a square matrix with the cross product of the relative abundances of taxa i,j. 
get_rel_ab <- function(vec, taxa){
  # make sure vectors are the same lenght
  stopifnot(length(vec) == length(taxa)) 
  # make a vector of relative abundances
  rel.ab <- vec / sum(vec) 
  # make an empty square matrix
  Nij <- matrix(0,
                length(vec),
                length(vec))
  # nested for loops to populate matrix i,j, with cross products of relative abundances of taxa i, and j, respectively.  
  for (i in 1:length(vec)){
    for (j in 1:length(vec)){
      Nij[i,j] <- rel.ab[i]*rel.ab[j]
    }
  }
  # add row and columns names from taxa vector
  dimnames(Nij) <- list(taxa, taxa)
  Nij
}

# function to scale continuous variable, x, to [min, max]
scalexy <- function(x,
                    min,
                    max){
  ((max - min) / (max(x) - min(x))) *
    (x - min(x)) + min
}

# function to set probabilities of forbidden taxa to 0
# matrix must have taxa for column names 
# f.taxa is a character vector of taxa names
rm_niche <- function(matrix,
                     f.taxa){
  for(name in 
      (colnames(matrix)[colnames(matrix) %in% f.taxa])){
    matrix[,name] <- 0
  }
  matrix
}

# function to get fw_measures and dominant eigenvalue
# requires functions from stability_fxns_sauve.R
# returns the stability metric, as well as standard food web measures
get_measures <- function(matr, # probability matrix
                         s2, # initial value, see stability() function from Sauve
                         trials = 1, # number of bernoulli trials to run
                         scale.Jij = FALSE, # scale interactions based on body size?
                         correlate.Jij = FALSE, # correlate positive and negative interactions?
                         correlate.value = 0.7){
  result <- list()
  for(i in 1:trials){
    A <- rm_cycle(b_trial(matr))
    J <- jacobian_binary(A)
    if(scale.Jij == TRUE){
      J = scale.Jij(J)
    }
    if(correlate.Jij == TRUE){
      J = correlate.Jij(J,
                        correlate.value = correlate.value)
    }
    stab = stability(J, s2 = s2)
    fw_meas <- Get.web.stats(A)
    result[[i]] <- c(stab = stab, fw_meas) 
  }
  return(result)
}

# function to scale interaction strengths by body size
# used internally in get_measures()
scale.Jij <- function(J){
  # J is a Jacobian matrix where elements Jij = interaction strength
  # e.g. object from jacobian_binary()
  element.rank <- rank(1:nrow(J)) / nrow(J)
  scale.rank <- scalexy(element.rank,
                        min = 0.25,
                        max = 1.25)
  scale.grid <- expand.grid(rev(scale.rank),
                            scale.rank)
  scale.vector <- scale.grid$Var1 * scale.grid$Var2
  scale.matrix <- matrix(scale.vector,
                         nrow = nrow(J),
                         ncol = ncol(J))
  scale.J <- scale.matrix * J
  return(scale.J)
}

# function to correlate positive and negative interactions by body size
# used internally in get_measures()
correlate.Jij <- function(J,
                          correlate.value = 0.7){
  # J is a JAcobian matrix where elements Jij = interaction strength
  # e.g. object from jacobian_binary()
  negative.index <- which(J < 0,
                          arr.ind = TRUE)
  positive.index <- cbind(negative.index[,2],
                          negative.index[,1])
  positive.strength <- abs(J[negative.index] * correlate.value)
    #runif(n = nrow(negative.index), min = 0.6, max = 0.8))
  # 0.7 comes from Montoya et al. 2009
  J[positive.index] <- positive.strength
  return(J)
}

# used internally in get_measures()
rm_cycle <- function(A, diag = TRUE){
  # if A eats B and B eats A, randomly remove one link
  # Default removes diagonal
  # to leave diagonal as is set diag = FALSE
  for(i in 1:nrow(A)){
  for(j in 1:ncol(A)){
    if(diag == FALSE){
      if(i == j){
        next
      }
    }
    if(A[i,j] == 1 & A[j,i] == 1){
      if(runif(1) < 0.5){
        A[i,j] <-  0
      }
      else{
        A[j,i] <- 0
      }
    }
  }
}
  return(A)
}

# food web functions modified from Owen Petchey
# original functions were found on thetrophiclink.org
# and are available on github here: https://github.com/opetchey/dumping_ground/blob/master/random_cascade_niche/FoodWebFunctions.r
Get.web.stats <- function(web, which.stats=1:4, real.web=NA){
  if(!is.list(web)){
    result <- c(S = length(web[,1]),
                L = sum(web),
                C = sum(web)/length(web[,1])^2)
    
    if(sum(which.stats==1)==1){
      BITUC <- Bottom.Intermediate.Top(web, proportion=T)$Proportions.of.each
      result <- c(result,
                  B = BITUC[1],
                  I = BITUC[2],
                  T = BITUC[3]
                  )
    }
  }
  
  if(is.list(web)){
    for(i in 1:length(web)){
      
      t.result <- c(S = length(web[[i]][,1]),
                    L = sum(web[[i]]),
                    C = sum(web[[i]])/length(web[[i]][,1])^2)
      
      if(sum(which.stats==1)==1){
        BITUC <- Bottom.Intermediate.Top(web[[i]], proportion=T)$Proportions.of.each
        t.result <- c(t.result,
                      B = BITUC[1],
                      I = BITUC[2],
                      T = BITUC[3]
                      )
      }
      
      if(i==1)
        result <- t.result
      if(i>1)
        result <- rbind(result, t.result)
    }
  }    
  result
}

# Cannibals <- function(web){
#   length(dimnames(web)[[1]][diag(web)==1]) / length(web[1,])
# }
## return the proportion of bottom, intermediate, top,
## unconnected and purely cannibalistic species
Bottom.Intermediate.Top <- function(web, proportion=TRUE, check.web=FALSE){
  
  ## find the names and numbers of BIT, unconnected and pure cannibals
  names.all <- 1:length(web[1,])
  dimnames(web) <- list(names.all, names.all)
  names.Bottom <-  names.all[apply(web, 2, sum)==0 & apply(web, 1, sum)!=0]#pulls out names where column sum ==0, and where row sum == >0. e.g. an alga does not consume anyhting (column sum == 0) but is consumed (row sum = >0)
  number.Bottom <- length(names.Bottom)
  names.Top <-   names.all[apply(web, 2, sum)!=0 & apply(web, 1, sum)==0]
  number.Top <- length(names.Top)
  names.Unconnected <- names.all[apply(web, 2, sum)==0 & apply(web, 1, sum)==0]
  number.Unconnected <- length(names.Unconnected)
  number.Pure.cannibals <- 0
  names.Pure.cannibals <- character(0)
  for(i in 1:length(web[,1]))
    if(web[i,i]!=0 && sum(web[-i,i]+web[i,-i])==0){
      if(number.Pure.cannibals==0){
        names.Pure.cannibals <- dimnames(web)[[1]][i]
        number.Pure.cannibals <- 1
      }
      else{
        names.Pure.cannibals <- c(names.Pure.cannibals, dimnames(web)[[1]][i])
        number.Pure.cannibals <- number.Pure.cannibals + 1
      }
    }
  names.Intermediate <- dimnames(web)[[1]][is.na(match(dimnames(web)[[1]],
                                                       c(names.Bottom, names.Top, names.Unconnected, names.Pure.cannibals)))]
  number.Intermediate <- length(web[1,])-number.Bottom-number.Top-number.Unconnected-number.Pure.cannibals
  
  if(proportion)
    result <- list(Bottom=names.Bottom,
                   Intermediate=names.Intermediate,
                   Top=names.Top,
                   Unconnected=names.Unconnected,
                   Pure.cannibals=names.Pure.cannibals,
                   Proportions.of.each=c(number.Bottom, number.Intermediate, number.Top, number.Unconnected, number.Pure.cannibals)/length(web[1,]))
  if(!proportion)
    result <- list(Bottom=names.Bottom,
                   Intermediate=names.Intermediate,
                   Top=names.Top,
                   Unconnected=names.Unconnected,
                   Pure.cannibals=names.Pure.cannibals,
                   Proportions.of.each=c(number.Bottom, number.Intermediate, number.Top, number.Unconnected, number.Pure.cannibals))
  result
}
