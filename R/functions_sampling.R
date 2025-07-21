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
                                             seed = NULL){
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
    
    random_uniform <- candidate_events %>%
      group_by(sender, replace) %>%
      summarize(r = runif(1), .groups = 'drop')
    candidate_events <- candidate_events %>%
      left_join(random_uniform, by = c("sender", "replace"))
    candidate_events <- candidate_events %>%
      mutate(waitingTimes = ifelse(sampled == 0, getWaitingTime(r, rate), NA))
    candidate_events <- candidate_events[,-6]
    
    minTime <- min(candidate_events$waitingTimes, na.rm = T)
    selectedEventID <- which(candidate_events$waitingTimes == minTime)
    
    if(length(selectedEventID)>1) {
      eventsTot.u <- arrange(rbind(creation.events.u,deletion.events.u), time)
      support.constrain.u = computeSupportConstrain(eventsTot.u, net1, actors.u$label)
      support.constrain.creation.u <- support.constrain.u[eventsTot.u$replace == 1]
      support.constrain.deletion.u <- support.constrain.u[eventsTot.u$replace == 0]
      
      # Need to check the choice probabilities 
      if(c(1) %in% candidate_events$replace[selectedEventID] ){
        creationChoice <- getDyNAMChoices(formulas_creation[[2]],
                                          initial_parameters_creation[[2]],
                                          panel_net, creation.events.u,
                                          support.constrain.creation.u,
                                          actors.u,rem_environment,
                                          rem_net,
                                          time1 = time1, time2 = time2)
      }else{ creationChoice = NULL}
      if(c(0) %in% candidate_events$replace[selectedEventID]){
        deletionChoice <- getDyNAMChoices(formulas_deletion[[2]], 
                                          initial_parameters_deletion[[2]],
                                          panel_net, deletion.events.u,
                                          support.constrain.deletion.u,
                                          actors.u,rem_environment,
                                          rem_net,
                                          time1 = time1, time2 = time2)
      }else{deteleChoice = NULL}
      
      
      auxDf <- candidate_events %>%
        mutate(senderId = sapply(sender, labelTorow,net1))
      auxDf <- auxDf %>%
        mutate(receiverId = sapply(receiver, labelTorow,net1))
      
      p <- apply(auxDf[selectedEventID,c("senderId","receiverId","replace")],1,
                 getChoiceProb, creationChoice, deletionChoice)
      
      selectedEventID <- sample(selectedEventID,1,prob = p)
    }
    # update the network
    chosenEvent <- candidate_events[selectedEventID,]
    
    if(chosenEvent$replace == 0){
      deletion.events.u <- deletion.events.u[-which(deletion.events.u$sender == chosenEvent$sender &
                                                      deletion.events.u$receiver == chosenEvent$receiver),]
    }else{
      creation.events.u <- creation.events.u[-which(creation.events.u$sender == chosenEvent$sender &
                                                      creation.events.u$receiver == chosenEvent$receiver),]
    }
    
    panel_net[ labelTorow(chosenEvent$sender,net1),labelTorow(chosenEvent$receiver,net1) ] <- chosenEvent$replace
    
    # update time
    simulationTime <- simulationTime + minTime
    
    # add event to the list
    simulatedEvents[nrow(simulatedEvents) + 1,] <- c(simulationTime, chosenEvent$sender, chosenEvent$receiver, chosenEvent$replace)
    
    # update candidate list
    candidate_events[selectedEventID, "sampled"] <- 1
    candidate_events <- candidate_events[,-6]
    
    nSampledEvents <- nSampledEvents + 1   
  }
  
  simulatedEvents$time <- as.numeric(simulatedEvents$time)
  simulationTime <- as.numeric(simulationTime)
  
  endTime <- simulationTime +
    mean(c(simulatedEvents$time, NA) -
           c(0, simulatedEvents$time), na.rm = T)
  
  simulatedEvents$time <- time1 + 
    (time2 - time1) * (simulatedEvents$time/endTime)
  
  return(simulatedEvents)
}
