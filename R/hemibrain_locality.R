##########################################################################################
################################ Overlap Locality Score ##################################
##########################################################################################

#' Get useful metrics for each neuron's axon and dendrite
#'
#' @description Calculates/retreives useful metrics for the axon and dendrites of a neuron, and its
#' split.
#'
#' @inheritParams hemibrain_overlap_locality
#' @param locality if TRUE \code{hemibrain_overlap_locality} is called and the
#' overlap locality score will be returned in addition to other metrics.
#'
#' @return a named vector of overlap_locality scores
#' @seealso \code{\link{hemibrain_overlap_locality}}
#' @export
hemibrain_compartment_metrics <- function(x, resample = 1000, delta = 1000, locality = TRUE, ...){
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
compartment_metrics <- function(x, resample = 1000, delta = 1000, locality = TRUE){

  # Axon-dendrite split?
  if(!(sum(x$d$Label==2)&sum(x$d$Label==3))){
    warning("Axon / dendrite missing")
  }

  # Synapses
  syns = tryCatch(hemibrain_extract_synapses(x), error = function(e) NULL)
  axon.pre = sum(syns$prepost==0&syns$Label==2)
  dend.pre = sum(syns$prepost==0&syns$Label==3)
  axon.post = sum(syns$prepost==1&syns$Label==2)
  dend.post = sum(syns$prepost==1&syns$Label==3)
  total.pre = sum(syns$prepost==0)
  total.post = sum(syns$prepost==1)

  # Segregation
  si = x$AD.segregation.index
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
  met = data.frame(
             total.pre = nullToNA(total.pre),
             total.post = nullToNA(total.post),
             axon.pre = nullToNA(axon.pre),
             dend.pre = nullToNA(dend.pre),
             axon.post = nullToNA(axon.post),
             dend.post = nullToNA(dend.post),
             axon.length= nullToNA(axon.length),
             dend.length= nullToNA(dend.length),
             pd.length= nullToNA(pd.length),
             segregation_index = nullToNA(si),
             overlap_locality = nullToNA(locality))
  met
}


#' Calculate the overlap score between a neurons axon dendrite
#'
#' @description Calculates an overlap score using \code{overlap_score} between a neuron's axon and dendrite,
#' if axon and dendrite have been marked out in the neuron object, e.g. by using \code{flow_centrality}.
#'
#' @inheritParams flow_centrality
#' @inheritParams overlap_score
#' @param resample stepsize to which to resample neurons. If set to \code{NULL},
#' neurons are not resampled.
#'
#' @return a named vector of overlap_locality scores
#' @seealso \code{\link{overlap_score}}
#' @export
hemibrain_overlap_locality <- function(x, resample = 1000, delta = 1000, ...){
  x = nat::as.neuronlist(x)
  l.scores = nat::nlapply(x, overlap_locality, resample = resample, delta = delta, ...)
  unlist(l.scores)
}

# hidden
overlap_locality <- function(x, resample = 1000, delta = 1000){

  # Axon-dendrite split?
  if(!(sum(x$d$Label==2)&sum(x$d$Label==3))){
    stop("Axon / dendrite missing")
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
  oa = overlap_score(axon.res, axon.res, delta = delta)
  od = overlap_score(dend.res, dend.res, delta = delta)
  omax = max(oa,od)
  onorm = oad/omax

  # Return
  onorm

}


#' Generate a connectivity matrix based on euclidean distance between points
#'
#' @description Generates an 'overlap matrix' of overlap scores between neurons in the \code{output.neurons} and \code{input.neurons} pools.
#' For every point in a given neuron in \code{output.neurons}, a distance score is calculated to every point in a neuron in \code{input.neurons}.
#' The sum of this score is added to the final output matrix. The score is calculated as \code{e(-d^2/2δ^2)}, where d is the euclidean distance between the two points,
#' and δ is the expected distance in um that is considered 'close'. It is recommended that the user resamples neurons before use, using \code{\link{resample}}.
#'
#' @param output.neurons first set of neurons
#' @param input.neurons second set of neurons
#' @param delta the distance at which a synapse might occur
#'
#' @examples
#' \dontrun{
#' # Calculate how much some neurons overlap with one another
#' ## Example requires the package nat.flybrains
#' Cell07PNs_overlap = overlap(output.neurons = nat.flybrains::Cell07PNs,
#' input.neurons = nat.flybrains::Cell07PNs)
#'
#' ## Plot the results
#' heatmap(Cell07PNs_overlap)
#' }
#' @return a matrix of overlap scores
#' @seealso \code{\link{hemibrain_overlap_locality}}
#' @export
overlap_score <- function(output.neurons, input.neurons, delta = 100){
  output.neurons = nat::as.neuronlist(output.neurons)
  input.neurons = nat::as.neuronlist(input.neurons)
  score.matrix = matrix(0,nrow = length(output.neurons), ncol = length(input.neurons))
  rownames(score.matrix) = names(output.neurons)
  colnames(score.matrix) = names(input.neurons)
  for (n in 1:length(output.neurons)){
    a = nat::xyzmatrix(output.neurons[[n]])
    input.neurons.d = nat::nlapply(input.neurons, nat::xyzmatrix, .progress = "none")
    s = sapply(input.neurons.d, function(x)sum(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2))))
    # Score similar to that in Schlegel et al. 2015
    score.matrix[n,] = s
  }
  score.matrix
}
