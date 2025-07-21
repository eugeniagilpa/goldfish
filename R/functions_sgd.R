##################### ###
#
# Goldfish package
# Functions to generate samples of sequences given the panel data.
# Differences for sampling for competition models or not are cosidered.
#
##################### ###


#' sgd
#' 
#' @param initialParameters
#' @param fixedParameters
#' @param nChains
#' @param formulas
#' @param formulasType
#' @param statistics
#' @param time1
#' @param time2
#' @param subModelTypes
#' @param environmentsChains
#' @param supportConstrain
#' @param dumping
#' @param eps
#' @param max.it
#' @param impSampling
#' @param chainSamplePRaw
#' @param nSizeBatch
#'  
#' @return new parameter estimators, standard errors and variances
#' @noRd
#'
sgd <- function(initialParameters,
                fixedParameters = NULL,
                nChains,
                formulas,
                formulasType,
                statistics,
                time1,
                time2,
                subModelTypes,
                environmentsChains,
                supportConstrain,
                dumping = 0.001,
                eps = 0.01,
                max.it = 100,
                impSampling = FALSE,
                chainSamplePRaw = 0,
                nSizeBatch = 1
                ){
  
  # get not fixed parameters
  fixed_positions <- lapply(fixedParameters, function(x) {
    if (is.null(x)) return(FALSE)
    pos <- which(!is.null(x) & !is.na(x))
    if (length(pos) == 0) return(FALSE)
    pos
  })
  
  for(i in 1:max.it){ #Iteration on the SGD algorithm
    # browser()
    
    # 1. Select the batch of samples that are going to be used in the update. # 
    # Sample one of the chains with importance sampling weights
    if(impSampling){
      sampleIndex <- sample(1:nChains, size = nSizeBatch, prob = chainSamplePRaw)
    }else{
      sampleIndex <- sample(1:nChains, size = nSizeBatch)
    }
    
    sampleIndex <- sort(sampleIndex)

    chainP <- chainSamplePRaw[sampleIndex]/sum(chainSamplePRaw[sampleIndex])
    # Compute score for that chain with the current values of initialParameters
    score <- vector(mode = "list", length = length(sampleIndex))
    infMatrixComplete <- vector(mode = "list", length = length(sampleIndex))
    
    calcProbInit <- lapply(1:length(initialParameters),
                           function(x) set_estimation_opt(
                             initial_parameters = initialParameters[[x]],
                             fixed_parameters = fixedParameters[[x]],
                             max_iterations = 0L,
                             engine = "default"
                           ))    
    # 2. Compute statistics for the batch sampled # 
    res <- lapply(seq_along(sampleIndex),  
                  function(i) Map(
                    function(formula, sub_model, preproc_init, op_list, control_est) {
                      estimate_dynam(
                        x = formula,
                        sub_model = sub_model,
                        preprocessing_init = preproc_init,
                        data = environmentsChains[[sampleIndex[i]]],
                        control_estimation = control_est,
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
                    statistics[[sampleIndex[i]]],
                    supportConstrain[[sampleIndex[i]]],
                    calcProbInit
                  ))  


    for( iChain in seq_along(sampleIndex)){
      score[[iChain]] <- lapply(res[[iChain]] , "[[", "finalScore")
      infMatrixComplete[[iChain]] <- lapply(res[[iChain]] , "[[", "finalInformationMatrix")
      
      if(impSampling){
        score[[iChain]] <- lapply(score[[iChain]], `*`,
                                  chainP[iChain])
        infMatrixComplete[[iChain]] <- lapply(infMatrixComplete[[iChain]], `*`,
                                              chainP[iChain])
        
      }else{
        score[[iChain]] <- lapply(score[[iChain]], `*`,
                                  1/nSizeBatch)
        infMatrixComplete[[iChain]] <- lapply(infMatrixComplete[[iChain]], `*`,
                                              1/nSizeBatch)
      }
    }
    if(!i%% 10)    cat("Iteration ", i,". Sample Index sgd routine: ", sampleIndex, "\n")
    
    
    totScore <- lapply(1:length(formulas), function(i) Reduce('+',lapply(score, "[[", i)) )
    step <-  lapply(totScore,`*`,dumping)
    
    
    
    step <- Map(function(vec, fix) {
      if (!identical(fix, FALSE)) vec[fix] <- 0
      vec
    }, step, fixed_positions)
    
    stop <-  norm(unlist(step),type = "2")
    if(!i%%10) cat("Iteration ", i,". Stop sgd routine: ", stop, "\n")
    # print(stop)
    if(stop<eps) {
      cat("stop < eps in sgd routine \n")
      break
    }
    initialParameters <- lapply(1:length(initialParameters),
                                function(x) initialParameters[[x]] + step[[x]] )
  }
  
  totInfMatrixComplete <- lapply(1:length(formulas),
                                 function(i) Reduce('+',lapply(infMatrixComplete, "[[", i)))
  onlyNonFixedInfMatrix <-  Map(function(mat, fix) {
    if (identical(fix, FALSE)) return(mat)
    mat[-fix, -fix, drop = FALSE]
  }, totInfMatrixComplete, fixed_positions)
  invTotInfMatrixComplete <- lapply(onlyNonFixedInfMatrix,solve)
  
  totScoreNonFixed <- Map(function(vec, fix) {
    if (identical(fix, FALSE)) return(vec)
    vec[-fix]
  }, totScore, fixed_positions)
  
  variances <- lapply(1:length(initialParameters),
                      function(i) invTotInfMatrixComplete[[i]]%*%
                        (totScoreNonFixed[[i]]%*%t(totScoreNonFixed[[i]])) %*% invTotInfMatrixComplete[[i]]) # following section 3.3 of Ruth (2024)
  # function(i) solve(totInfMatrixComplete[[i]] - totScore[[i]]%*%t(totScore[[i]]))) #information matrix ruth 2.2
  
  
  variances <- Map(add_na_rowcol, variances, fixed_positions)
  stdErrors <- Map(std_error_fixed_params, variances, fixed_positions)

  return(list("initialParameters" = initialParameters, 
              "stdErrors" = stdErrors,
              "variances" = variances))
}


