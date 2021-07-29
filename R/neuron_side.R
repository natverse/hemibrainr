# Functions for calculating brain side

#' Calculate which side of the brain a neuron is on
#'
#' @description These functions make use of the symmetrised brainspace \code{JRC2018F} to determine the 'side' of the brain on which a neuron's root node, typically the soma, falls, i.e. \code{soma_side}.
#' You can also calculate a 'synapse side index', defined as:  (number of synapses on the right - number on the left) / total synapses.
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param brain wthe brainspace what \code{x} is in.
#' @param bound the dividing line that separates left and right in e \code{JRC2018F}.
#'
#' @return a \code{nat::neuronlist} object with either 'side' or synapse_side_index updated in its metadata.
#'
#' @examples
#' \donttest{
#' library(hemibrainr)
#' library(nat.jrcbrains)
#'
#' # Read a neuron from the hemibrain
#' neurons = neuprint_read_neurons(1796818119)
#'
#' # Split into axon and dendrite
#' neurons.flow = flow_centrality(neurons)
#'
#' # Calculate soma side
#' neurons.flow = soma_side(neurons.flow, brain = "JRCFIB2018F")
#'
#' # Calculate synapse side index
#' neurons.flow = synapse_side_index(neurons, brain = "JRCFIB2018F")
#'
#' # See side data
#' neurons.flow[,]
#'
#'}
#' @export
#' @seealso \code{\link{hemibrain_reroot}}
#' @export
soma_side <-function(x,
                     brain = "JRC2018F",
                     bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2) UseMethod("soma_side")

#' @export
soma_side.neuronlist <- function(x,
                                 brain = "JRC2018F",
                                 bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2){
  if(brain!="JRC2018F"){
    y = java_xform_brain(x, reference = "JRC2018F", sample = brain, .parallel = FALSE, verbose = FALSE, OmitFailures = FALSE, progress.rjava=TRUE)
  }else{
    y = x
  }
  leftsomas = unlist(nat::nlapply(y, soma_side.neuron, brain = "JRC2018F", bound = bound))
  leftsomas = leftsomas[!is.na(leftsomas)]
  x[names(leftsomas),"side"] = leftsomas
  x
}

#' @export
soma_side.neuron <- function(x,
                      brain = "JRC2018F",
                      bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2){
  if(brain!="JRC2018F"){
    x = java_xform_brain(x, reference = "JRC2018F", sample = brain, .parallel = FALSE, verbose = FALSE, OmitFailures = FALSE, progress.rjava=TRUE)
  }
  r = nat::rootpoints(x)
  if(is.numeric(r)){
    position = nat::xyzmatrix(x$d[r,])
    ifelse(position[,"X"]>bound,"left","right")
  }else{
    NA
  }
}

#' @export
synapse_side_index <-function(x,
                              brain = "JRC2018F",
                     bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2) UseMethod("synapse_side_index")

#' @export
synapse_side_index.neuronlist <- function(x,
                                 brain = "JRC2018F",
                                 bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]-nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2){
  if(brain!="JRC2018F"){
    y = java_xform_brain(x, reference = "JRC2018F", sample = brain, .parallel = FALSE, verbose = FALSE, OmitFailures = FALSE, progress.rjava=TRUE)
  }else{
    y = x
  }
  ssi = nat::nlapply(y, synapse_side_index.neuron, brain = "JRC2018F", bound = bound)
  ssi = do.call(rbind,ssi)
  rownames(ssi) = names(y)
  x[rownames(ssi),colnames(ssi)] = ssi
  x
}

#' @export
synapse_side_index.neuron <- function(x,
                                      brain = "JRC2018F",
                               bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2){
  if(brain!="JRC2018F"){
    x = java_xform_brain(x, reference = "JRC2018F", sample = brain, .parallel = FALSE, verbose = FALSE, OmitFailures = FALSE, progress.rjava=TRUE)
  }
  x$connectors$side = ifelse(nat::xyzmatrix(x$connectors)[,"X"]>bound,"left","right")
  agg = aggregate(list(count = x$connectors$connector_id),
                  list(side = x$connectors$side,
                       prepost = x$connectors$prepost),
                  function(x) length(unique(x)))
  pre.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==0)$count)
  pre.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==0)$count)
  presynapse_side_index = (pre.right-pre.left)/(pre.right+pre.left)
  post.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==1)$count)
  post.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==1)$count)
  postsynapse_side_index = (post.right-post.left)/(post.right+post.left)
  df = data.frame(postsynapse_side_index=postsynapse_side_index, presynapse_side_index=presynapse_side_index)
  if(!is.null(x$connectors$Label)){
    x$connectors$Label = standard_compartments(x$connectors$Label)
    agg = aggregate(list(count = x$connectors$connector_id),
                    list(side = x$connectors$side,
                         prepost = x$connectors$prepost,
                         Label = x$connectors$Label),
                    function(x) length(unique(x)))
    # axons
    pre.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==0 & agg$Label == "axon")$count)
    pre.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==0 & agg$Label == "axon")$count)
    axon_presynapse_side_index = (pre.right-pre.left)/(pre.right+pre.left)
    post.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==0 & agg$Label == "axon")$count)
    post.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==0 & agg$Label == "axon")$count)
    axon_postsynapse_side_index = (post.right-post.left)/(post.right+post.left)
    # dendrite
    pre.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==0 & agg$Label == "dendrite")$count)
    pre.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==0 & agg$Label == "dendrite")$count)
    dendrite_presynapse_side_index = (pre.right-pre.left)/(pre.right+pre.left)
    post.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==1 & agg$Label == "dendrite")$count)
    post.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==1 & agg$Label == "dendrite")$count)
    dendrite_postsynapse_side_index = (post.right-post.left)/(post.right+post.left)
    df = data.frame(postsynapse_side_index=postsynapse_side_index, presynapse_side_index=presynapse_side_index,
                    axon_postsynapse_side_index=axon_postsynapse_side_index, axon_presynapse_side_index=axon_presynapse_side_index,
                    dendrite_postsynapse_side_index=dendrite_postsynapse_side_index, dendrite_presynapse_side_index=dendrite_presynapse_side_index)
  }
}

# hidden
emptytozero <- function(v){
  if(!length(v)){
    0
  }else{
    v
  }
}
