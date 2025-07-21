##################### ###
#
# Goldfish package
# Some utility functions useful in
# various parts of the code, specific
# to the emDyNAM algorithm
#
##################### ###


#' unfixed_params
#'
#' @param params list
#' @param fixedParameters list
#'
#' @return Subset of parameters not fixed
#' @noRd
#'
#' @examples
#' \donttest{
#' params <- list(c(0,0),c(0,0))
#' fixedParameters <- list(NULL,c(NA,1))
#' unfixed_params(params,fixedParameters)
#' }
unfixed_params <- function(params,fixedParameters){
  fixed_positions <- lapply(fixedParameters, function(x) {
    if (is.null(x)) return(FALSE)
    pos <- which(!is.null(x) & !is.na(x))
    if (length(pos) == 0) return(FALSE)
    pos
  })
  
  params <- Map(function(params, fix) {
    if (!identical(fix, FALSE)) params <- params[-fix]
    params
  }, params, fixed_positions)
  return(params)
}

#' make_fix_effects
#' 
#' Given a (list of) formula(s), looks up for a certain string (effect) and
#' creates (a list of) fixed paramters following `golfish` standards. All 
#' parameters are left unfixed, except for the positions of `char`, that are
#' fixed to `val`.  
#'
#' @param fmla (list of) formula(s) 
#' @param char string 
#' @param val double
#'
#' @return Fixed parameters formated for goldfish estimation
#' @noRd
#'
#' @examples
#' \donttest{
#' fmla <- list(net ~ 1 + outdeg + indeg)
#' make_fix_effects(fmla,"fixedParameters,"out",0)
#' }
make_fix_effects <- function(fmla, char, val) {
  rhs <- as.character(fmla)[3]
  terms <- trimws(unlist(strsplit(rhs, "\\+")))
  idx <- grep(char, terms)
  if (length(idx) == 0) {
    return(NULL)
  } else {
    result <- rep(NA, length(terms))
    result[idx] <- val
    return(result)
  }
}


#' clipping
#' 
#' Given a vector, truncates values larger than the `quan` quantile.
#'
#' @param x vector of doubles
#' @param quant quantile ($\in[0,1]$) 
#'
#' @return vector of doubles truncated
#' @noRd
#'
#' @examples
#' \donttest{
#' clipping(seq(0,1,0.1),0.5)
#' }
clipping <- function(x,quant = 0.9){
  thr <- quantile(x , quant)
  x[x > thr] <- thr
  return(x)
}



#' get_weights
#' 
#' Compute importance sampling weigths
#'
#' @param proposal vector of log-proposals
#' @param loglik dataframe of log-likelihoods
#' @param f function transformation of weigths
#'
#' @return vector of normalized importance sampling weights
#' @noRd
#'
get_weights <- function(proposal, loglik, f = function(x) x ,...){
  auxRelWeight <-rowSums(loglik) -min(rowSums(loglik)) - proposal +min(proposal)
  w <- exp(auxRelWeight)
  w <- f(w,...)
  w <- w / sum(w) 
  return(w)
}



#' stratified_resample
#' 
#' Stratified resampling method
#' 
#' @references Singh, R., Mangat, N.S. (1996). Stratified Sampling. In: Elements of Survey Sampling. Kluwer Texts in the Mathematical Sciences, vol 15. Springer, Dordrecht. https://doi.org/10.1007/978-94-017-1404-4_5
#' @references J. D. Hol, T. B. Schon and F. Gustafsson, "On Resampling Algorithms for Particle Filters," 2006 IEEE Nonlinear Statistical Signal Processing Workshop, Cambridge, UK, 2006, pp. 79-82, doi: 10.1109/NSSPW.2006.4378824.
#'
#' @param cum_weights vector of cumulative weights
#' 
#' @return vector with resampled indices
#' @noRd
#'
stratified_resample <- function(cum_weights) {
  N <- length(cum_weights)
  u <- (runif(N) + 0:(N-1)) / N # Stratified uniform samples
  indices <- findInterval(u, cum_weights) + 1
  return(indices)
}


#' residual_resample
#' 
#' Residual resampling method
#'
#' @references J. D. Hol, T. B. Schon and F. Gustafsson, "On Resampling Algorithms for Particle Filters," 2006 IEEE Nonlinear Statistical Signal Processing Workshop, Cambridge, UK, 2006, pp. 79-82, doi: 10.1109/NSSPW.2006.4378824. 
#'
#' @param weights vector of weights
#' 
#' @return vector with resampled indices
#' @noRd
#'
residual_resample <- function(weights) {
  N <- length(weights)
  weights <- weights / sum(weights)  # Ensure weights are normalized
  
  # Step 1: Deterministic part
  num_copies <- floor(N * weights)
  residual <- weights * N - num_copies
  indices <- rep(1:N, num_copies)
  
  # Step 2: Stochastic part
  R <- N - length(indices)  # Number of remaining particles to sample
  if (R > 0) {
    residual <- residual / sum(residual)  # Normalize residuals
    cumulative <- cumsum(residual)
    u <- runif(R)
    extra_indices <- findInterval(u, cumulative) + 1
    indices <- c(indices, extra_indices)
  }
  
  return(indices)
}




#' labelTorow
#' 
#' Utility function to transform vector of node labels to the rows in the given network
#'
#' @param vec vector node labels
#' @param net adjacency matrix of network of interest
#' 
#' @return newVec with row number corresponding to the nodes
#' @noRd
#'
labelTorow <- function(vec,net){
  # vec: vector with labels that we want to transform to row/col number
  newVec = sapply(1:length(vec), function(i) which(as.integer(colnames(net)) == vec[i]) )
  return(newVec)
}


#' getChoiceProb
#' 
#' Utility function to get the choice probability for a given an event
#' TODO: restructure to be more generic and not crea-del focused
#'
#' @param row row of a sequence of events: replace == 0 implies deletion of events
#' @param creationChoice matrix of creation choice probabilities
#' @param deletionChoice matrix of deletion choice probabilities
#' 
#' @return choice probability of the event
#' @noRd
#'
getChoiceProb <- function(row, creationChoice, deletionChoice) {
  if (row[3] == 0) {
    return(deletionChoice[row[1], row[2]])
  } else {
    return(creationChoice[row[1], row[2]])
  }
}



#' getWaitingTime
#' 
#' Utility function to get a waiting time (exponential distribution)
#'
#' @param r random uniform number in [0,1]
#' @param lambda parameter of the exponential distribution
#' 
#' @return waiting time to next event
#' @noRd
#'
getWaitingTime <- function(r, lambda){
  - log(1-r) / lambda 
}



#' getWaitingTime
#' 
#' Utility function to get the supportConstrain (opportunityList) for a sequence of events
#' TODO: restructure to be more generic and not crea-del focused
#'
#' @param eventsTot data.frame sequence of events
#' @param net1 matrix initial network
#' @param label node labels of the network
#' 
#' @return list with support constraints for creation and deletion events
#' @noRd
#'
computeSupportConstrain <- function(eventsTot, net1, label){
  
  supportConstrain <- vector(mode = "list", length = nrow(eventsTot))
  auxNet <- net1
  
  for(i in 1:nrow(eventsTot)){
    sender <- eventsTot$sender[i]
    receiver <- eventsTot$receiver[i]
    
    if(eventsTot$replace[i] == 1){
      isAvailableChoice <- which(auxNet[labelTorow(sender,net1),]==0 & label != sender)
    }else{
      isAvailableChoice <- which(auxNet[labelTorow(sender,net1),]==1)
    }
    supportConstrain[[i]] <- isAvailableChoice
    auxNet[labelTorow(sender,net1),labelTorow(receiver,net1)] <- eventsTot$replace[i]
  }
  supportConstrainCrea <- supportConstrain[eventsTot$replace == 1]
  supportConstrainDel <- supportConstrain[eventsTot$replace == 0]
  
  return(list("creation" = supportConstrainCrea, "deletion" = supportConstrainDel))
  
}


