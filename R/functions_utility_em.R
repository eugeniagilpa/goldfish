##################### ###
#
# Goldfish package
# Some utility functions useful in
# various parts of the code, specific
# to the emDyNAM algorithm
#
##################### ###


## Random utility functioins #-----


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
unfixed_params <- function(params,fixedParameters) {
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
clipping <- function(x,quant = 0.9) {
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
get_weights <- function(proposal, loglik, f = function(x) x ,...) {
  auxRelWeight <- rowSums(loglik) - min(rowSums(loglik)) - proposal + min(proposal)
  w <- exp(auxRelWeight)
  w <- f(w,...)
  w <- w / sum(w)
  return(w)
}



#' stratified_resample
#'
#' Stratified resampling method
#'
#' @references Singh, R., Mangat, N.S. (1996). Stratified Sampling
#' In: Elements of Survey Sampling. Kluwer Texts in the Mathematical
#' Sciences, vol 15. Springer, Dordrecht.
#' https://doi.org/10.1007/978-94-017-1404-4_5
#' @references J. D. Hol, T. B. Schon and F. Gustafsson, "On Resampling
#' Algorithms for Particle Filters," 2006 IEEE Nonlinear Statistical Signal
#' Processing Workshop, Cambridge, UK, 2006, pp. 79-82,
#' doi: 10.1109/NSSPW.2006.4378824.
#'
#' @param cum_weights vector of cumulative weights
#'
#' @return vector with resampled indices
#' @noRd
#'
stratified_resample <- function(cum_weights) {
  N <- length(cum_weights)
  u <- (runif(N) + 0:(N - 1)) / N # Stratified uniform samples
  indices <- findInterval(u, cum_weights) + 1
  return(indices)
}


#' residual_resample
#'
#' Residual resampling method
#'
#' @references J. D. Hol, T. B. Schon and F. Gustafsson,
#'  "On Resampling Algorithms for Particle Filters," 2006 IEEE
#' Nonlinear Statistical Signal Processing Workshop, Cambridge,
#' UK, 2006, pp. 79-82, doi: 10.1109/NSSPW.2006.4378824
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



#' add_na_rowcol
#'
#' @param mat matrix
#' @param fixed_positions vector of positions where to add NAs
#'
#' @return matrix with NA in the fixed positions (row and column)
#' @noRd
#'
add_na_rowcol <- function(mat, fixed_positions) {
  if (identical(fixed_positions, FALSE) || length(fixed_positions) == 0) return(mat)

  n <- nrow(mat)
  new_n <- n + length(fixed_positions)
  new_mat <- matrix(NA, nrow = new_n, ncol = new_n)
  all_idx <- seq_len(new_n)
  insert_idx <- fixed_positions
  orig_idx <- setdiff(all_idx, insert_idx)
  new_mat[orig_idx, orig_idx] <- mat
  return(new_mat)
}




#' std_error_fixed_params
#'
#' @param mat matrix
#' @param fixed_positions vector of positions where to add NAs
#'
#' @return vector of std errors with 0 in the fixed positions
#' @noRd
#'
std_error_fixed_params <- function(mat, fixed_positions) {
  d <- sqrt(diag(mat))
  if (!is.null(fixed_positions) && length(fixed_positions) > 0) {
    d[fixed_positions] <- 0
  }
  d
}

## Probability-related utility functions #----

#' getChoiceProb
#'
#' Utility function to get the choice probability for a given an event
#' TODO: restructure to be more generic and not crea-del focused
#'
#' @param row row of a sequence of events: replace == 0 implies
#' deletion of events
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
getWaitingTime <- function(r, lambda) {
  - log(1 - r) / lambda
}

#test

#' computeSupportConstrain
#'
#' Utility function to get the supportConstrain (opportunityList)
#' for a sequence of events
#' TODO: restructure to be more generic and not crea-del focused
#'
#' @param eventsTot data.frame sequence of events
#' @param net1 matrix initial network
#' @param label node labels of the network
#'
#' @return list with support constraints for creation and deletion events
#' @noRd
#'
computeSupportConstrain <- function(eventsTot, net1, label) {

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

#test

#' make_environment_from_list
#'
#' Former loadDataFast
#'
#' @param panel_events events from the imputation in the panel data (data.frame)
#' @param rem_events data.frame relational events
#' @param panel_data logical, if TRUE, the panel data goldfish object
#' is created, else only actors in the environment
#' @param rem_data logical, if TRUE, the fb.net goldfish object
#' is created, else, only actors in the environment
#' @param net matrix initial panel network
#' @param time1 start time
#' @param time2 end time
#' @param actors data.frame actors, must have column present
#' @param presentUpdate data.frame with updates of the presence 
#' of actors panel_data
#'
#' @return goldfish make_data object with panel and rem data for estimation
#' @noRd
#'
make_environment_from_list <- function(panel_events = NULL,
                                       rem_events = NULL,
                                       panel_data = TRUE,
                                       rem_data = TRUE,
                                       net = NULL,
                                       time1 = NULL,
                                       time2 = NULL,
                                       actors = NULL,
                                       presentUpdate = NULL
                                       ) {
  if(panel_data){
    if(!is.character(panel_events$sender)) panel_events$sender=as.character(panel_events$sender)
    if(!is.character(panel_events$receiver)) panel_events$sender=as.character(panel_events$receiver)
  }

    datObj <- getData(actors = actors,
                      presentUpdate = presentUpdate,
                      panel_events = panel_events,
                      rem_events = rem_events,
                      net = net,
                      time1 = time1,
                      time2 = time2,
                      panel_data = panel_data,
                      rem_data = rem_data
                      )

  return(datObj)
}



#' getData
#' 
#' @param panel_events events from the imputation in the panel data (data.frame)
#' @param rem_events data.frame relational events
#' @param panel_data logical, if TRUE, the panel data goldfish object is created
#' else only actors in the environment
#' @param rem_data logical, if TRUE, the fb.net goldfish object is created,
#' else, only actors in the environment
#' @param net matrix initial panel network
#' @param time1 start time
#' @param time2 end time
#' @param actors data.frame actors, must have column present
#' @param presentUpdate data.frame with updates of the presence of actors
#' rem_data
#'  
#' @return goldfish make_data object with panel and rem data for estimation,
#' events from sequence "i"
#' @noRd
#'
getData <- function(panel_events = NULL,
                    rem_events = NULL,
                    panel_data = TRUE,
                    rem_data = TRUE,
                    panel_net = NULL,
                    time1 = NULL,
                    time2 = NULL,
                    actors = NULL,
                    presentUpdate = NULL) {

  actors <- make_nodes(actors)
  actors_rem <- make_nodes(actors[, c("label", "present")])
  actors_rem$present <- FALSE
  if(!is.null(presentUpdate)) actors_rem <- link_events(actors_rem, presentUpdate, attribute = "present")
  actors$present <- TRUE

  if(panel_data){
    panel_net <- make_network(net, nodes = actors, directed = TRUE)
    panel_net <- link_events(panel_net, change_events = panel_events, nodes = actors)
    panel_dependent <- make_dependent_events(events = panel_events,
                                             nodes = actors,
                                             default_network = panel_net)
  }
  if(rem_data){
    rem_events <- rem_events[rem_events$time <= time2, ]
    rem_net <- make_network(nodes = actors_rem, directed = T)
    rem_net <- link_events(rem_net, change_events = rem_events, nodes = actors_rem)
    rem_dependent <- make_dependent_events(events = rem_events,
                                           nodes = actors_rem,
                                           default_network = rem_net)
    rem_dependent <- rem_dependent[rem_dependent$time > time1, ]
  }

  if (panel_data) {
    if (rem_data) {
      data_dynam <- make_data(
        actors = actors,
        actors_rem = actors_rem,
        presentUpdate = presentUpdate,
        panel_net = panel_net,
        panel_dependent = panel_dependent,
        panel_events = panel_events,
        rem_net = rem_net,
        rem_dependent = rem_dependent,
        rem_events = rem_events
      )
    }else {
      data_dynam <- make_data(
        actors = actors,
        actors_rem = actors_rem,
        presentUpdate = presentUpdate,
        panel_net = panel_net,
        panel_dependent = panel_dependent,
        panel_events = panel_events
        )
    }
  }else {
    if (rem_data) {
      data_dynam <- make_data(
        actors = actors,
        actors_rem = actors_rem,
        presentUpdate = presentUpdate,
        rem_net = rem_net,
        rem_dependent = rem_dependent,
        rem_events = rem_events
      )
    }else {
      data_dynam <- make_data(
        actors = actors,
        actors_rem = actors_rem,
        presentUpdate = presentUpdate
      )
    }
  }
  return(data_dynam)
}




#' getDyNAMRates
#'
#' @param rate_formula formula rate model
#' @param rate_params vector of parameters for the rate model
#' @param net.u network
#' @param events.u data.frame sequence of events
#' @param actors.u actors
#' @param rem_environment environment relational events
#' @param time1 start time
#' @param time2 end time
#'
#' @return node poisson rates
#' @noRd
#'
getDyNAMRates <- function(rate_formula,
                          rate_params,
                          net.u,
                          events.u = NULL,
                          actors.u,
                          rem_environment = NULL,
                          time1 = 0,
                          time2 = 1
                          ) {


  if (nrow(events.u) == 0 || is.null(events.u)) {
    events.u <- data.frame(time = (time1 + (time2 - time1) / 2),
                           sender = actors.u$label[1],
                           receiver = actors.u$label[2],
                           replace = Inf)
  }


  panel_net <- make_network(net.u, nodes = actors.u, directed = TRUE)
  panel_net <- link_events(panel_net, change_events = events.u, nodes = actors.u)


  dep_events <- make_dependent_events(
    events.u,
    nodes = actors.u, default_network = panel_net)
  rate_formula[[2]] <- as.name("dep_events")


  data <- make_data(
    actors.u = actors.u,
    friend.net = panel_net,
    depEvents = dep_events
  )

  if (!is.null(rem_environment)) {
    for (name in ls(envir = rem_environment, all.names = TRUE)) {
      assign(name, get(name, envir = rem_environment), envir = data)
    }
  }


  rate_stats <- estimate_dynam(
    rate_formula,
    sub_model = "rate",
    data = data,
    preprocessing_only = TRUE,
    verbose = FALSE,
    progress = FALSE
  )

  # transform matrices to vectors
  node_stats <- apply(rate_stats$initialStats, 3, function(x) apply(x, 1, max))

  # calculate node poisson rates
  node_rates <- as.vector(exp(rate_params[1] +
                                rowSums(t(t(node_stats) * rate_params[-1]))
  )
  )

  return(node_rates)
}






#' getDyNAMChoices
#'
#' @param rate_formula formula rate model
#' @param rate_params vector of parameters for the rate model
#' @param net.u network
#' @param events.u data.frame sequence of events
#' @param actors.u actors
#' @param rem_environment environment relational events
#' @param time1 start time
#' @param time2 end time
#' @param support.constrain.u list of support constraints for the events
#'
#' @return node poisson rates
#' @noRd
#'
getDyNAMChoices <- function(choice_formula,
                            choice_params,
                            net.u,
                            events.u = NULL,
                            actors.u,
                            rem_environment = NULL,
                            time1 = 0,
                            time2 = 1,
                            support.constrain.u = NULL
                            ) {

  if (nrow(events.u) == 0 || is.null(events.u)) {
    events.u <- data.frame(time = (time1 + (time2 - time1) / 2),
                           sender = actors.u$label[1],
                           receiver = actors.u$label[2],
                           replace = Inf)
  }


  panel_net <- make_network(net.u, nodes = actors.u, directed = TRUE)
  panel_net <- link_events(panel_net, change_events = events.u, nodes = actors.u)


  dep_events <- make_dependent_events(
    events.u,
    nodes = actors.u,
    default_network = panel_net)
  choice_formula[[2]] <- as.name("dep_events")


  data <- make_data(
    actors.u = actors.u,
    panel_net = panel_net,
    dep_events = dep_events
  )

  if (!is.null(rem_environment)) {
    for (name in ls(envir = rem_environment, all.names = TRUE)) {
      assign(name, get(name, envir = rem_environment), envir = data)
    }
  }


  choice_stats <- estimate_dynam(choice_formula,
                                sub_model = "choice",
                                data = data,
                                preprocessing_only = TRUE,
                                verbose = FALSE,
                                progress = FALSE
  )

  nEffects <- dim(choice_stats$initialStats)[3]
  weightedStats <- array(0, dim = c(nrow(actors.u), nrow(actors.u), nEffects))
  for (i in 1:nEffects){
    weightedStats[,,i] <- choice_stats$initialStats[,,i] * choice_params[i]
  }
  objFun <- exp(apply(weightedStats, 1:2, sum))
  diag(objFun) <- 0
  choiceProbabilities <- objFun * (rowSums(objFun)^(-1))

  tieChoices <- choiceProbabilities

  return(tieChoices)
}





#' mleMC
#'
#' @param indexCore
#' @param splitIndicesPerCore
#' @param formulas list of formulas
#' @param subModelTypes list of sub-model types
#' @param statistics list of preprocessed goldfish objects
#' @param supportConstrain list of support constraints for each formula
#' @param time1 start time
#' @param time2 end time
#' @param sampleIndexUnique vector of resampled sequences
#' @param mleNeeded vector of indices to compute MLE
#' @param data list of goldfish environments
#' @param fixed_parameters list of fixed parameters for the estimation
#'
#' @return node poisson rates
#' @noRd
#'
mleMC  <- function(indexCore,
                   splitIndicesPerCore,
                   formulas = formulas,
                   subModelTypes = subModelTypes,
                   statistics = statistics,
                   supportConstrain = supportConstrain,
                   time1 = time1,
                   time2 = time2,
                   sampleIndexUnique = sampleIndexUnique,
                   mleNeeded = mleNeeded,
                   data = environmentsChains,
                   fixed_parameters = NULL
                   ) {

  indicesCore <- splitIndicesPerCore[[indexCore]]
  resMLE <- vector("list", length(indicesCore))
  for (aux_sampleID in seq_along(indicesCore)) {
    sampleID <- indicesCore[aux_sampleID]
    id <- sampleIndexUnique[mleNeeded[sampleID]]

    res <- Map(
      function(formula, sub_model, preproc_init, op_list, fix) {
        estimate_dynam(
          x = formula,
          sub_model = sub_model,
          preprocessing_init = preproc_init,
          data = environmentsChains[[id]],
          control_estimation = set_estimation_opt(engine = "default",
                                                  fixed_parameters = fix),
          control_preprocessing = set_preprocessing_opt(
            start_time = time1,
            end_time = time2,
            opportunities_list = op_list
          ),
          verbose = FALSE,
          progress = FALSE
        )
      },
      formulas,
      subModelTypes,
      statistics[[id]],
      supportConstrain[[id]],
      fixedParameters
    )

    resMLE[[aux_sampleID]] <- res
  }

  return(resMLE)
}
