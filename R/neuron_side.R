# Functions for calculating brain side

#' Calculate which side of the brain a neuron is on
#'
#' @description These functions make use of the symmetrised brainspace \code{JRC2018F} to determine the 'side' of the brain on which a neuron's root node, typically the soma, falls, i.e. \code{soma_side}.
#' You can also calculate a 'synapse side index', defined as:  (number of synapses on the right - number on the left) / total synapses.
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param brain wthe brainspace what \code{x} is in.
#' @param bound the dividing line that separates left and right in e \code{JRC2018F}.
#' @param transform whether to transform \code{x} (when \code{TRUE}) when the given brain is not \code{JRC2018F}
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
#' @seealso \code{\link{hemibrain_reroot}}
#' @export
#' @rdname soma_side
soma_side <-function(x,
                     brain = nat.flybrains::JRC2018F,
                     bound = NULL,
                     transform = FALSE) UseMethod("soma_side")

#' @export
#' @rdname soma_side
soma_side.neuronlist <- function(x,
                                 brain = nat.flybrains::JRC2018F,
                                 bound = NULL,
                                 transform = FALSE){

  ss = sort_side(x=x, brain=brain, bound=bound, transform=transform)
  y = ss$y
  bound = ss$bound
  leftsomas = unlist(nat::nlapply(y, soma_side.neuron, brain = "JRC2018F", bound = bound))
  leftsomas = leftsomas[!is.na(leftsomas)]
  x[names(leftsomas),"side"] = leftsomas
  x
}

#' @export
#' @rdname soma_side
soma_side.neuron <- function(x,
                             brain = nat.flybrains::JRC2018F,
                             bound = NULL,
                             transform = FALSE){
  if(is.null(bound)){
    ss = sort_side(x=x, brain=brain, bound=bound, transform=transform)
    y = ss$y
    bound = ss$bound
  }else{
    y=x
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
#' @rdname soma_side
synapse_side_index <-function(x,
                              brain = nat.flybrains::JRC2018F,
                              bound = NULL,
                              transform = FALSE) UseMethod("synapse_side_index")

#' @export
#' @rdname soma_side
synapse_side_index.neuronlist <- function(x,
                                          brain = nat.flybrains::JRC2018F,
                                          bound = NULL,
                                          transform = FALSE){
  ss = sort_side(x=x, brain=brain, bound=bound, transform=transform)
  y = ss$y
  bound = ss$bound
  # soma side
  leftsomas = unlist(nat::nlapply(y, soma_side.neuron, bound = bound))
  leftsomas = leftsomas[!is.na(leftsomas)]
  x[names(leftsomas),"side"] = leftsomas
  # synapse side
  ssi = nat::nlapply(y, synapse_side_index.neuron, bound = bound)
  ssi = do.call(rbind,ssi)
  rownames(ssi) = names(y)
  x[rownames(ssi),colnames(ssi)] = ssi
  # return
  x
}

#' @export
#' @rdname soma_side
synapse_side_index.neuron <- function(x,
                                      brain = nat.flybrains::JRC2018F,
                                      bound = NULL,
                                      transform = FALSE){
  if(is.null(bound)){
    ss = sort_side(x=x, brain=brain, bound=bound, transform=transform)
    y = ss$y
    bound = ss$bound
  }else{
    y=x
  }
  x$connectors$side = ifelse(nat::xyzmatrix(x$connectors)[,"X"]>bound,"left","right")
  if(is.nrowlength(x$connectors)){
    agg = stats::aggregate(list(count = x$connectors$connector_id),
                    list(side = x$connectors$side,
                         prepost = x$connectors$prepost),
                    function(x) length(unique(x)))
  }else{
    agg = data.frame(count = 0, side = c("right","left"), prepost = c(0,1,1,0))
  }
  pre.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==0)$count)
  pre.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==0)$count)
  presynapse_side_index = (pre.right-pre.left)/(pre.right+pre.left)
  post.right = emptytozero(subset(agg, agg$side=="right" & agg$prepost==1)$count)
  post.left = emptytozero(subset(agg, agg$side=="left" & agg$prepost==1)$count)
  postsynapse_side_index = (post.right-post.left)/(post.right+post.left)
  df = data.frame(postsynapse_side_index=postsynapse_side_index, presynapse_side_index=presynapse_side_index)
  if(!is.null(x$connectors$Label)){
    x$connectors$Label = standard_compartments(x$connectors$Label)
    if(length(x$connectors$Label)){
      agg = stats::aggregate(list(count = x$connectors$connector_id),
                      list(side = x$connectors$side,
                           prepost = x$connectors$prepost,
                           Label = x$connectors$Label),
                      function(x) length(unique(x)))
    }else{
      agg1 = data.frame(count = 0, side = c("right","left"), prepost = c(0,1,1,0), Label = c("axon"))
      agg2 = data.frame(count = 0, side = c("right","left"), prepost = c(0,1,1,0), Label = c("dendrite"))
      agg = rbind(agg1,agg2)
    }
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

# hidden
sort_side <- function(x, brain, bound = NULL, transform = FALSE){
  boundnull = is.null(bound)
  if(is.character(brain)){
    sample = brain
  }else{
    sample = brain$name
  }
  if(sample!="JRC2018F"){
    if(transform){
      y = nat.templatebrains::xform_brain(x, reference = "JRC2018F", sample = sample, .parallel = TRUE, verbose = FALSE, OmitFailures = FALSE)
      if(boundnull) bound = (nat::boundingbox(nat.flybrains::JRC2018F.surf)[2,1]+nat::boundingbox(nat.flybrains::JRC2018F.surf)[1,1])/2
    }else{
      y = x
      if(boundnull) bound = (nat::boundingbox(brain)[2,1]+nat::boundingbox(brain)[1,1])/2
    }
  }else{
    y = x
    if(boundnull) bound = (brain$BoundingBox[2,1]+brain$BoundingBox[1,1])/2
  }
  list(y = y, brain = brain, sample = sample, bound = bound, transform = transform)
}
