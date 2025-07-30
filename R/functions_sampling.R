##################### ###
#
# Goldfish package
# Functions to generate samples of sequences given the panel data.
# Differences for sampling for competition models or not are cosidered.
#
##################### ###

### Get chain from model functions #------

#' get_chain_from_competition_model
#' 
#' @param formulas_creation
#' @param formulas_deletion
#' @param initial_parameters_creation 
#' @param initial_parameters_deletion 
#' @param candidate_events total candidate events (maybe to remove?)
#' @param creation.events.u creation events
#' @param deletion.events.u deletion events
#' @param time1 start time
#' @param time2 end time
#' @param actors.u actors
#' @param rem_environment relational data goldfish environment
#' @param net1 initial network
#' @param verbose logical, if TRUE prints additional information
#' @param seed integer, seed for random number generation
#'  
#' @return node poisson rates
#' @noRd
#'
get_chain_from_competition_model <- function(formulas_creation, 
                                             formulas_deletion,
                                             initial_parameters_creation,
                                             initial_parameters_deletion,
                                             candidate_events,
                                             creation.events.u,
                                             deletion.events.u,
                                             time1, 
                                             time2,
                                             actors.u,
                                             rem_environment,
                                             net1,
                                             verbose = F,
                                             seed = NULL,
                                             condProbs = T){
  if(!is.null(seed)) set.seed(seed)
  
  simulatedEvents <- data.frame(time = numeric(0), sender = numeric(0), receiver = numeric(0), replace = numeric(0))
  candidate_events$sampled  <- 0
  simulationTime <- 0
  nSampledEvents <- 0
  panel_net <- net1
  
  # loop through candidate events 
  while(any(candidate_events$sampled == 0)){
    
    nRemainingEvents <- sum(candidate_events$sampled == 0)
    candidate_events$rate <- NA

    # calculate rates and waiting times
    creationRates <- getDyNAMRates(formulas_creation[[1]],
                                   initial_parameters_creation[[1]],
                                   panel_net, creation.events.u, actors.u,
                                   rem_environment,
                                   rem_net,
                                   time1 = time1, time2 = time2)
    deletionRates <- getDyNAMRates(formulas_deletion[[1]],
                                   initial_parameters_deletion[[1]],
                                   panel_net, deletion.events.u, actors.u,
                                   rem_environment,
                                   rem_net,
                                   time1 = time1, time2 = time2)
    
    
    # allows that there are multiple creations between the same two actors
    creationIndexes <- subset(candidate_events, sampled == 0 & replace == 1)[, c("sender", "receiver")]
    deletionIndexes <- subset(candidate_events, sampled == 0 & replace == 0)[, c("sender", "receiver")]
    
    if(nrow(creationIndexes)>0){
      candidate_events[candidate_events$sampled == 0 & candidate_events$replace == 1,]$rate <-
        creationRates[labelTorow(creationIndexes$sender,net1)]
    }
    
    if(nrow(deletionIndexes)>0){    
      candidate_events[candidate_events$sampled == 0 & candidate_events$replace == 0,]$rate <- 
        deletionRates[labelTorow(deletionIndexes$sender,net1)]
    }
    
    aux <- unique(candidate_events[candidate_events$sampled == 0,c('sender','rate')])
    if(condProbs) {
      aux$rateP <- aux$rate / sum(aux$rate)
      selectedEventIDaux <- sample(1:nrow(aux),1,prob =  aux$rateP)
      proposalProb <- proposalProb  + log(aux$rateP[selectedEventIDaux])
    }else{
      selectedEventIDaux <- sample(1:nrow(aux),1,prob =  aux$rate)
      proposalProb <- proposalProb + log(aux$rate[selectedEventIDaux])
    }
    
    selectedEventID <- which(candidate_events$sender == aux[selectedEventIDaux,]$sender)
    selectedEventID <- candidate_events[selectedEventID,]
    selectedEventID <- selectedEventID[ selectedEventID$sampled==0,]
    selectedEventID <- selectedEventID[abs(selectedEventID$rate-aux[selectedEventIDaux,]$rate)<1e-16,]
    
    if(nrow(selectedEventID)>1) {
      # Need to check the choice probabilities 
      if(c(1) %in% selectedEventID$replace ){
        creationChoice <- getDyNAMChoices(choice_formula = formulas_creation[[2]],
                                          choice_params = initial_parameters_creation[[2]],
                                          net.u = net.u, 
                                          events.u = candidate_events, 
                                          actors.u = actors.u,,
                                          rem_environment = rem_environment,
                                          time1 = time1, 
                                          time2 = time2)
        
        
      }else{ creationChoice = NULL}
      if(c(0) %in% selectedEventID$replace){
        deletionChoice <- getDyNAMChoices(choice_formula = formulas_deletion[[2]],
                                          choice_params = initial_parameters_deletion[[2]],
                                          net.u = net.u, 
                                          events.u = candidate_events, 
                                          actors.u = actors.u,,
                                          rem_environment = rem_environment,
                                          time1 = time1, 
                                          time2 = time2)
      }else{deteleChoice = NULL}
      
      auxDf <- selectedEventID %>%
        mutate(senderId = sapply(sender, labelTorow,net1))
      auxDf <- auxDf %>%
        mutate(receiverId = sapply(receiver, labelTorow,net1))
      
      p <- apply(auxDf[,c("senderId","receiverId","replace")],1,
                 getChoiceProb, creationChoice, deletionChoice)
      if(condProbs) p<- p/sum(p)
      selectedEventIDaux <- sample(1:nrow(selectedEventID),1,prob = p)
      
      proposalProb <- proposalProb + log(p[selectedEventIDaux])
      
    }
    # update the network
    if(nrow(selectedEventID)>1){
      chosenEvent <- selectedEventID[selectedEventIDaux,]
    }else{
      chosenEvent <- selectedEventID
    }
    
    net.u[ labelTorow(chosenEvent$sender,net1),labelTorow(chosenEvent$receiver,net1) ] <- chosenEvent$replace

    simulationTime <- simulationTime + getWaitingTime(runif(1), chosenEvent$rate)
    
    # add event to the list
    simulatedEvents[nrow(simulatedEvents) + 1,] <- c(simulationTime, chosenEvent$sender, chosenEvent$receiver, chosenEvent$replace)
    
    # update candidate list
    matching_row <- which(apply(candidateEvents, 1, function(row) all(row == chosenEvent)))
    candidateEvents[matching_row, "sampled"] <- 1
    candidateEvents <- candidateEvents[,-6]
    
    nSampledEvents <- nSampledEvents + 1   
    
  }
  
  
  simulatedEvents$time <- as.numeric(simulatedEvents$time)
  simulationTime <- as.numeric(simulationTime)
  
  
  endTime <- simulationTime +
    mean(c(simulatedEvents$time, NA) -
           c(0, simulatedEvents$time), na.rm = T)
  
  
  simulatedEvents$time <- time1 + 
    (time2 - time1) * (simulatedEvents$time/endTime)
  
  returnList <- list(simulatedEvents = simulatedEvents,
                     proposalProb = proposalProb)
  return(returnList)
}


#' get_chain_from_competition_model_mc
#' 
#' Multi-core auxiliary function for get_chain_from_competition_model
#' 
#' @param indexCore
#' @param splitIndicesPerCore
#' @param formulas_creation
#' @param formulas_deletion
#' @param initial_parameters_creation 
#' @param initial_parameters_deletion 
#' @param candidate_events total candidate events (maybe to remove?)
#' @param creation.events.u creation events
#' @param deletion.events.u deletion events
#' @param time1 start time
#' @param time2 end time
#' @param actors.u actors
#' @param rem_environment relational data goldfish environment
#' @param net1 initial network
#' @param verbose logical, if TRUE prints additional information
#' @param seed integer, seed for random number generation
#' @param condProbs logical, if TRUE computes conditional probabilities
#'  
#' @return node poisson rates
#' @noRd
#'
get_chain_from_competition_model_mc <- function(indexCore,
                                       splitIndicesPerCore,
                                       formulas_creation, # 1 rate, 2 choice
                                       formulas_deletion,
                                       initial_parameters_creation,
                                       initial_parameters_deletion,
                                       candidate_events,
                                       creation.events.u,
                                       deletion.events.u,
                                       time1, 
                                       time2,
                                       actors.u,
                                       rem_environment,
                                       net1,
                                       verbose = F,
                                       seed=NULL,
                                       condProbs = T
){
  print(splitIndicesPerCore)
  print(indexCore)
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resChain <- vector("list", length(indicesCore))
  resProposalProb <- vector("list", length(indicesCore))
  

  res <- lapply(seed[indicesCore],get_chain_from_competition_model_mc,
                formulas_creation, 
                formulas_deletion,
                initial_parameters_creation,
                initial_parameters_deletion,
                candidate_events,
                creation.events.u,
                deletion.events.u,
                time1, 
                time2,
                actors.u,
                rem_environment,
                net1,
                verbose,
                seed,
                condProbs)
  
  resChain <- lapply(res, function(x) x$simulatedEvents)
  resProposalProb <- lapply(res, function(x) x$proposalProb)
  resultT <- list(resChain= resChain,
                  resProposalProb = resProposalProb)
  return(resultT)
  
}

#' get_chain_competition
#' 
#' @param net1
#' @param net2
#' @param time1 start time
#' @param time2 end time
#' @param labels kabeks of actors
#' @param seed integer, seed for random number generation
#'  
#' @return sequence between net1 and net2
#' @noRd
#'
get_chain_competition <- function(net1, 
                                  net2,
                                  time1, 
                                  time2,
                                  labels,
                                  seed = NULL){
  
  if(!is.null(seed)) set.seed(seed)
  
  edgesC <-  which(net2 > net1, arr.ind = T)
  edgesD <-  which(net2 < net1, arr.ind = T)
  
  # draw relative times from two constant poisson parameters
  # this equals a random uniform sampling
  if(!is.null(seed)) set.seed(seed)
  timesC <- runif(nrow(edgesC))
  timesD <- runif(nrow(edgesD))
  
  edges <- data.frame(rbind(edgesC,edgesD))
  names(edges) <- c("sender","receiver") 
  edges$replace <- c(rep(1,nrow(edgesC)),rep(0,nrow(edgesD)))
  edges$time <- c(timesC, timesD)
  
  edges <- edges[sample(1:nrow(edges)),]  
  
  simulationTime <- sum(edges$time)
  
  endTimecreation <- simulationTime + mean(timesC)
  endTimedeletion <- simulationTime + mean(timesD)
  
  for(j in 2:nrow(edges)) edges$time[j] <- edges$time[j]+edges$time[j-1]
  
  edges$time[edges$replace == 1]  <- time1 +
    (time2 - time1) * (edges$time[edges$replace == 1]/endTimecreation)
  edges$time[edges$replace == 0] <- time1 + 
    (time2 - time1) * (edges$time[edges$replace == 0]/endTimedeletion)
  
  edges <- data.frame("time" = edges$time,
                      "sender" = labels[edges$sender],
                      "receiver" = labels[edges$receiver],
                      "replace" = edges$replace)
  
  return(edges)  
}





