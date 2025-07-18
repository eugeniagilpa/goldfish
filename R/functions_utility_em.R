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


