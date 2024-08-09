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
#' @examples
#' \donttest{
#'
#' # Choose neurons
#' ## In this case some antennal lobe local neurons
#' al.local.neurons = c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#' "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#' "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#' "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#' "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#' "1794037618", "5813018729", "2013333009")
#'
#' # use a smaller group in examples
#' al.local.neurons=al.local.neurons[1:5]
#'
#' # Get neurons
#' neurons = neuprint_read_neurons(al.local.neurons)
#'
#' # Convert to microns
#' neurons = hemibrainr::scale_neurons(neurons)
#'
#' # Split
#' neurons = hemibrain_flow_centrality(neurons)
#'
#' # Calculate overlap and other metrics
#' hcm = hemibrain_compartment_metrics(neurons,
#' resample = NULL,
#' just.leaves = TRUE, # just look how much the leaves of skeletons overlap, quicker (default)
#' delta = 5)
#'
#'}
#' @return a \code{data.frame} of metrics broken down by putative axon and dendrite, and their overlap_locality scores
#' @seealso \code{\link{hemibrain_overlap_locality}}
#' @export
hemibrain_compartment_metrics <- function(x, resample = 10, delta = 62.5, locality = FALSE, ...){
  mets = nat::nlapply(x, compartment_metrics, resample = resample, delta = delta, locality = locality, ...)
  mets = mets[!unlist(sapply(mets, is.null))]
  mets.df = as.data.frame(do.call(rbind, mets))
  mets.df$id = names(mets)
  if(nrow(mets.df)==length(x)){
    unik = setdiff(colnames(x[,]), colnames(mets.df))
    cbind(x[,unik],mets.df)
  }else{
    warning("Some neurons were dropped ...")
    mets.df
  }
}

# hidden
compartment_metrics <- function(x, resample = 10, delta = 62.5, locality = FALSE, ...){
  syns = x$connectors # tryCatch(hemibrain_extract_synapses(x), error = function(e) NULL)
  labd = ifelse("label" %in% colnames(x$d), "label", "Label")
  lab = ifelse("label" %in% colnames(syns), "label", "Label")
  both = (sum(x$d[[labd]]%in%c(2,"axon"))&sum(x$d[[labd]]%in%c(3,"dendrite")))
  # Axon-dendrite split?
  if(!sum(x$d[[labd]]%in%c(2,"axon"))&!sum(x$d[[labd]]%in%c(3,"dendrite"))){
    warning("axon/dendrite missing")
    total_length = tryCatch(summary(x)$cable.length, error = function(e) NA)
    total_outputs = NA
    total_inputs = NA
    axon_outputs = NA
    dend_outputs = NA
    pd_outputs = NA
    pn_outputs = NA
    axon_inputs = NA
    dend_inputs = NA
    pd_inputs = NA
    pn_inputs = NA
    total_outputs_density = NA
    total_inputs_density = NA
    axon_outputs_density = NA
    dend_outputs_density = NA
    axon_inputs_density = NA
    dend_inputs_density = NA
    axon_length = NA
    dend_length = NA
    pd_length = NA
    pn_length = NA
    axon_width = NA
    dend_width = NA
    pd_width = NA
    pn_width = NA
    si = NA
    locality.score = NA
  }else{

    # Synapses
    axon_outputs = tryCatch(sum(syns$prepost==0&syns[[lab]]%in%c(2,"axon")), error = function(e) 0)
    dend_outputs = tryCatch(sum(syns$prepost==0&syns[[lab]]%in%c(3,"dendrite")), error = function(e) 0)
    axon_inputs = tryCatch(sum(syns$prepost==1&syns[[lab]]%in%c(2,"axon")), error = function(e) 0)
    dend_inputs = tryCatch(sum(syns$prepost==1&syns[[lab]]%in%c(3,"dendrite")), error = function(e) 0)
    pd_outputs = tryCatch(sum(syns$prepost==0&syns[[lab]]%in%c(4,"primary_dendrite")), error = function(e) 0)
    pn_outputs = tryCatch(sum(syns$prepost==0&syns[[lab]]%in%c(7,"primary_neurite")), error = function(e) 0)
    pd_inputs = tryCatch(sum(syns$prepost==1&syns[[lab]]%in%c(4,"primary_dendrite")), error = function(e) 0)
    pn_inputs = tryCatch(sum(syns$prepost==1&syns[[lab]]%in%c(7,"primary_neurite")), error = function(e) 0)
    total_outputs = tryCatch(sum(syns$prepost==0), error = function(e) 0)
    total_inputs = tryCatch(sum(syns$prepost==1), error = function(e) 0)

    # Segregation
    if(both){
      si = round(tryCatch(x$AD.segregation.index, error = function(e) NA),6)
      if(locality){
        locality.score = tryCatch(overlap_locality(x, resample = resample, delta = delta, ...), error = function(e) NA)
      }else{
        locality.score = NA
      }
    }else{
      locality.score = NA
      si = NA
    }

    # Cable length
    axon_length = tryCatch(summary(axonic_cable(x))$cable.length, error = function(e) 0)
    dend_length = tryCatch(summary(dendritic_cable(x))$cable.length, error = function(e) 0)
    pd_length = tryCatch(summary(primary_dendrite_cable(x))$cable.length, error = function(e) 0)
    pn_length = tryCatch(summary(primary_neurite_cable(x))$cable.length, error = function(e) 0)
    total_length = tryCatch(summary(x)$cable.length, error = function(e) NA)

    # Cable width
    axon_width = tryCatch(cable_width(axonic_cable(x)), error = function(e) 0)
    dend_width = tryCatch(cable_width(dendritic_cable(x)), error = function(e) 0)
    pd_width = tryCatch(cable_width(primary_dendrite_cable(x)), error = function(e) 0)
    pn_width = tryCatch(cable_width(primary_neurite_cable(x)), error = function(e) 0)
  }

  # Assemble
  met = tryCatch(data.frame(
             total_outputs = nullToNA(total_outputs),
             total_inputs = nullToNA(total_inputs),
             axon_outputs = nullToNA(axon_outputs),
             dend_outputs = nullToNA(dend_outputs),
             axon_inputs = nullToNA(axon_inputs),
             dend_inputs = nullToNA(dend_inputs),
             pd_outputs = nullToNA(pd_outputs),
             pn_outputs = nullToNA(pn_outputs),
             pd_inputs = nullToNA(pd_inputs),
             pn_inputs = nullToNA(pn_inputs),
             total_outputs_density = nullToNA(total_outputs)/nullToNA(total_length),
             total_inputs_density = nullToNA(total_inputs)/nullToNA(total_length),
             axon_outputs_density = nullToNA(axon_outputs)/nullToNA(axon_length),
             dend_outputs_density = nullToNA(dend_outputs)/nullToNA(dend_length),
             axon_inputs_density = nullToNA(axon_inputs)/nullToNA(axon_length),
             dend_inputs_density = nullToNA(dend_inputs)/nullToNA(dend_length),
             total_length = nullToNA(total_length),
             axon_length= nullToNA(axon_length),
             dend_length= nullToNA(dend_length),
             pd_length= nullToNA(pd_length),
             pn_length=nullToNA(pn_length),
             axon_width= nullToNA(axon_width),
             dend_width= nullToNA(dend_width),
             pd_width= nullToNA(pd_width),
             pn_width=nullToNA(pn_width),
             segregation_index = nullToNA(si),
             overlap_locality = nullToNA(locality.score), stringsAsFactors = FALSE),
             error = function(e) NULL)
  if(!locality){
    met$overlap_locality = NULL
  }
  met = t(apply(met, 1, function(row) ifelse(is.nan(row),0,row)))
  met = apply(met,2,function(c) signif(c, digits = 6))
  t(as.data.frame(met, stringsAsFactors = FALSE))
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
                             delta = 62.5,
                             ...){

  # Axon-dendrite split?
  if(!(sum(x$d$Label%in%c(2,"axon"))&sum(x$d$Label%in%c(3,"dendrite")))){
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

  # Remove sub-trees with no synapses
  dend.res = prune_synapseless_branches(dend.res)
  axon.res = prune_synapseless_branches(axon.res)

  # Calculate compartment overlap scores
  oad = overlap_score_delta(axon.res, dend.res, delta = delta, ...)
  oda = overlap_score_delta(dend.res, axon.res, delta = delta, ...)
  oa = overlap_score_delta(axon.res, axon.res, delta = delta, ...)
  od = overlap_score_delta(dend.res, dend.res, delta = delta, ...)
  onorm = sum(oad,oda)/sum(oa,od)

  # Return
  onorm

}

# hidden, similar function now in nat
overlap_score_delta <- function(output.neurons, input.neurons, delta = 62.5, just.leaves = TRUE, max = exp(-delta^2/(2*delta^2)), normalise = FALSE){
  output.neurons = nat::as.neuronlist(output.neurons)
  input.neurons = nat::as.neuronlist(input.neurons)
  score.matrix = matrix(0,nrow = length(output.neurons), ncol = length(input.neurons))
  rownames(score.matrix) = names(output.neurons)
  colnames(score.matrix) = names(input.neurons)
  if(just.leaves){
    input.neurons.d = nat::nlapply(input.neurons, function(x) nat::xyzmatrix(x)[nat::endpoints(x),], .progress = "none")
  }else{
    input.neurons.d = nat::nlapply(input.neurons, nat::xyzmatrix, .progress = "none")
  }
  for (n in 1:length(output.neurons)){
    if(just.leaves){
      a = nat::xyzmatrix(output.neurons[[n]])[nat::endpoints(output.neurons[[n]]),]
    }else{
      a = nat::xyzmatrix(output.neurons[[n]])
    }
    if(normalise){
      s = sapply(input.neurons.d, function(x)
        lengthnorm(maxout(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)),max=max))
        )
    }else{
      s = sapply(input.neurons.d, function(x)
        sum(maxout(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)),max=max))
      ) # Score similar to that in Schlegel et al. 2015
    }
    score.matrix[n,] = s
  }
  score.matrix
}

#
projection_score.neuronlist <- function(x,
                                        resample = NULL,
                                        ...){
  y = nat::nlapply(x,projection_score.neuron, ...)
  data.frame(id = names(y), projection_score = unlist(y))
}

projection_score.neuron <- function(x,
                                    resample = NULL,
                                    ...){

  # Axon-dendrite split?
  if(!(sum(x$d$Label%in%c(2,"axon"))&sum(x$d$Label%in%c(3,"dendrite")))){
    warning("Axon / dendrite missing")
    return(NA)
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

  # Calculate midpoints
  axon.mid = colMeans(nat::xyzmatrix(axon.res))
  dend.mid = colMeans(nat::xyzmatrix(dend.res))

  # Euclidean distance between them
  euc.dist = stats::dist(rbind(axon.mid, dend.mid), method = "euclidean")[1]
  euc.dist

}




segregation_index.neuronlist <- function(x, ...){
  nat::nlapply(x, segregation_index.neuron, ...)
}

segregation_index.neuron <- function(x, ...){
  si <- segregation_index(x)
  x$AD.segregation.index <- si
  x
}

segregation_index <- function(neuron){
  nodes = neuron$d
  if(is.null(nodes)){
    return(NA)
  }else if (!nrow(nodes)){
    return(NA)
  }
  connectors = neuron$connectors
  if(is.null(connectors)){
    return(NA)
  }else if (!nrow(connectors)){
    return(NA)
  }
  label1 <- ifelse("Label"%in%colnames(nodes),"Label","label")
  label2 <- ifelse("Label"%in%colnames(connectors),"Label","label")
  dendrites = subset(nodes, nodes[[label1]] == 3)
  dendrites.post = nrow(subset(connectors, connectors$prepost == 1 & connectors[[label2]] %in% c("dendrite","3")))
  dendrites.pre = nrow(subset(connectors, connectors$prepost == 0  & connectors[[label2]] %in% c("dendrite","3")))
  dendrites.both = dendrites.post + dendrites.pre
  dendrites.pi = dendrites.post/dendrites.both
  dendrites.si = -(dendrites.pi * log(dendrites.pi) + (1 -dendrites.pi) * log(1 - dendrites.pi))
  if (is.nan(dendrites.si)) {
    dendrites.si = 0
  }
  axon = subset(nodes, nodes[[label1]] == 2)
  axon.post = nrow(subset(connectors, connectors$prepost == 1 & connectors[[label2]] %in% c("axon","2")))
  axon.pre = nrow(subset(connectors, connectors$prepost == 0  & connectors[[label2]] %in% c("axon","2")))
  axon.both = axon.post + axon.pre
  axon.pi = axon.post/axon.both
  axon.si = -(axon.pi * log(axon.pi) + (1 - axon.pi) * log(1 -axon.pi))
  if (is.nan(axon.si)) {
    axon.si = 0
  }
  entropy.score = (1/(dendrites.both + axon.both)) * ((axon.si *axon.both) + (dendrites.si * dendrites.both))
  both.comps = (dendrites.post + axon.post)/(dendrites.both +axon.both)
  control.score = -(both.comps * log(both.comps) + (1 - both.comps) *log(1 - both.comps))
  segregation.index = 1 - (entropy.score/control.score)
  if(is.na(segregation.index)||is.nan(segregation.index)) {
    segregation.index = 0
  }
  segregation.index
}








