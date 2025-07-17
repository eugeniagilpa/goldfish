### Auxiliar functions # ----
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
