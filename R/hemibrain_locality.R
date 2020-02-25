##########################################################################################
################################ Overlap Locality Score ##################################
##########################################################################################

#' Get useful metrics for each neuron's axon and dendrite
#'
#' @description Calculates/retrieves useful metrics for the axon and dendrites
#'   of a neuron, and its split.
#'
#' @inheritParams hemibrain_overlap_locality
#' @param locality if TRUE \code{\link{hemibrain_overlap_locality}} is called
#'   and the overlap locality score will be returned in addition to other
#'   metrics.
#'
#' @return a named vector of overlap_locality scores
#' @seealso \code{\link{hemibrain_overlap_locality}}
#' @export
hemibrain_compartment_metrics <- function(x, resample = 10, delta = 62.5, locality = TRUE, ...){
  mets = nat::nlapply(x, compartment_metrics, resample = resample, delta = delta, locality = locality, ...)
  mets.df = do.call(rbind, mets)
  if(nrow(mets.df)==length(x)){
    cbind(x[,],mets.df)
  }else{
    warning("Some neurons were dropped ...")
    mets.df
  }
}

# hidden
compartment_metrics <- function(x, resample = 10, delta = 62.5, locality = TRUE){

  # Axon-dendrite split?
  if(!(sum(x$d$Label==2)&sum(x$d$Label==3))){
    warning("Axon / dendrite missing")
  }

  # Synapses
  syns = tryCatch(hemibrain_extract_synapses(x), error = function(e) NULL)
  axon.outputs = tryCatch(sum(syns$prepost==0&syns$Label==2), error = function(e) NA)
  dend.outputs = tryCatch(sum(syns$prepost==0&syns$Label==3), error = function(e) NA)
  axon.inputs = tryCatch(sum(syns$prepost==1&syns$Label==2), error = function(e) NA)
  dend.inputs = tryCatch(sum(syns$prepost==1&syns$Label==3), error = function(e) NA)
  total.outputs = tryCatch(sum(syns$prepost==0), error = function(e) NA)
  total.inputs = tryCatch(sum(syns$prepost==1), error = function(e) NA)

  # Segregation
  si = tryCatch(x$AD.segregation.index, error = function(e) NA)
  if(locality){
    locality = tryCatch(overlap_locality(x, resample = resample, delta = delta), error = function(e) NA)
  }else{
    locality = NA
  }

  # Cable length
  axon.length = tryCatch(summary(axonic_cable(x))$cable.length, error = function(e) NA)
  dend.length = tryCatch(summary(dendritic_cable(x))$cable.length, error = function(e) NA)
  pd.length = tryCatch(summary(primary_dendrite_cable(x))$cable.length, error = function(e) NA)
  total.length = tryCatch(summary(x)$cable.length, error = function(e) NA)

  # Assemble
  met = tryCatch(data.frame(
             total.outputs = nullToNA(total.outputs),
             total.inputs = nullToNA(total.inputs),
             axon.outputs = nullToNA(axon.outputs),
             dend.outputs = nullToNA(dend.outputs),
             axon.inputs = nullToNA(axon.inputs),
             dend.inputs = nullToNA(dend.inputs),
             total.outputs.density = nullToNA(total.outputs)/nullToNA(total.length),
             total.inputs.density = nullToNA(total.inputs)/nullToNA(total.length),
             axon.outputs.density = nullToNA(axon.outputs)/nullToNA(axon.length),
             dend.outputs.density = nullToNA(dend.outputs)/nullToNA(dend.length),
             axon.inputs.density = nullToNA(axon.inputs)/nullToNA(axon.length),
             dend.inputs.density = nullToNA(dend.inputs)/nullToNA(dend.length),
             total.length = nullToNA(total.length),
             axon.length= nullToNA(axon.length),
             dend.length= nullToNA(dend.length),
             pd.length= nullToNA(pd.length),
             segregation_index = nullToNA(si),
             overlap_locality = nullToNA(locality)),
             error = function(e) NULL)
  met
}


#' Calculate the overlap score between a neurons axon dendrite
#'
#' @description Calculates an overlap score using an overlap score between a neuron's axon and dendrite,
#' if axon and dendrite have been marked out in the neuron object, e.g. by using \code{flow_centrality}.
#'
#' @inheritParams flow_centrality
#' @param delta the distance at which a synapse might occur
#' @param resample stepsize to which to resample neurons. If set to \code{NULL},
#' neurons are not resampled.
#'
#' @return a named vector of overlap_locality scores
#' @seealso \code{\link{hemibrain_compartment_metrics}}
#' @export
hemibrain_overlap_locality <- function(x, resample = 10, delta = 62.5, ...){
  x = nat::as.neuronlist(x)
  l.scores = nat::nlapply(x, overlap_locality, resample = resample, delta = delta, ...)
  unlist(l.scores)
}

# hidden
overlap_locality <- function(x,
                             resample = 10,
                             delta = 62.5){

  # Axon-dendrite split?
  if(!(sum(x$d$Label==2)&sum(x$d$Label==3))){
    warning("Axon / dendrite missing")
  }

  # Extract compartment cable
  axon = axonic_cable(x)
  dend = dendritic_cable(x)

  # Resample compartments
  if(!is.null(resample)){
    axon.res = nat::resample(axon, stepsize = resample)
    dend.res = nat::resample(dend, stepsize = resample)
  }else{
    axon.res = axon
    dend.res = dend
  }

  # Calculate compartment overlap scores
  oad = overlap_score(axon.res, dend.res, delta = delta)
  oda = overlap_score(dend.res, axon.res, delta = delta)
  oa = overlap_score(axon.res, axon.res, delta = delta)
  od = overlap_score(dend.res, dend.res, delta = delta)
  omax = max(oa,od)
  onorm = sum(oad,oda)/sum(oa,od)

  # Return
  onorm

}

# hidden, similar function now in nat
overlap_score <- function(output.neurons, input.neurons, delta = 62.5){
  output.neurons = nat::as.neuronlist(output.neurons)
  input.neurons = nat::as.neuronlist(input.neurons)
  score.matrix = matrix(0,nrow = length(output.neurons), ncol = length(input.neurons))
  rownames(score.matrix) = names(output.neurons)
  colnames(score.matrix) = names(input.neurons)
  for (n in 1:length(output.neurons)){
    a = nat::xyzmatrix(output.neurons[[n]])
    input.neurons.d = nat::nlapply(input.neurons, nat::xyzmatrix, .progress = "none")
    s = sapply(input.neurons.d, function(x)sum(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)))) # Score similar to that in Schlegel et al. 2015
    score.matrix[n,] = s
  }
  score.matrix
}
