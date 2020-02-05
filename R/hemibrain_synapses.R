############################################################################
################################ Synapses ##################################
############################################################################

#' Extract synapse location and labels from a neuron/neuronlist
#'
#' @description Extract a single data frame describing synapse types, partners,locations and
#' position on a neuron's axon/dendrite.
#' @inheritParams flow_centrality
#' @param prepost whether to get presynapses, postsynapses or both
#'
#' @return a \code{data.frame}
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Choose some known troble makers
#' id = c("452677169")
#'
#' # Read in these neurons
#' neuron = neuprintr::neuprint_read_neurons(id)
#'
#' # Re-root
#' neuron.flow = flow_centrality(neuron, polypre = TRUE,
#' mode = "centrifugal",
#' split = "distance")
#'
#' # Let's check that this worked
#' syns = hemibrain_extract_synapses(neuron.flow)
#'
#' # See result
#' plot3d_split(neuron.flow)
#' points3d(subset(syns,prepost==1), col = "cyan")
#' points3d(subset(syns,prepost==0), col = "red")
#'
#' }}
#' @export
#' @seealso \code{\link{flow_centrality}}
hemibrain_extract_synapses <- function(x,
                                       prepost = c("BOTH","PRE","POST"),
                                       ...){
  prepost = match.arg(prepost)
  if(nat::is.neuronlist(x)){
    syns = nat::nlapply(x,extract_synapses, ...)
    syns = do.call(rbind,syns)
  }else if (nat::is.neuron(x)){
    syns = extract_synapses(x)
  }else{
    stop("x must be a neuron or neuronlist object")
  }
  if(prepost=="PRE"){
    syns = syns[syns$prepost==0,]
  }else if (prepost=="POST"){
    syns = syns[syns$prepost==1,]
  }
  rownames(syns) = 1:nrow(syns)
  syns
}

# hidden
extract_synapses <-function(x){
  syn = x$connectors
  if(!nrow(syn)){
    warning("Neuron ", x$bodyid," has no synapses")
  }else{
    syn$bodyid = x$bodyid
    syn$skid = x$skid
    syn$Label = x$d$Label[match(syn$treenode_id,x$d$PointNo)]
  }
  syn
}
