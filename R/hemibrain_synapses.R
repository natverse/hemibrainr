############################################################################
################################ Synapses ##################################
############################################################################

#' Extract synapse location and labels from a neuron/neuronlist
#'
#' @description Extract a single data frame describing synapse/connection types,
#'   partners, locations and position on a neuron's axon/dendrite. You can
#'   either get all synapses returned or all unitary connections to a neuron's
#'   partners returned. Broken down by axon/dendrite (\code{Label}), and
#'   pre/postsynapses or pre/postsynaptic partners.
#' @inheritParams flow_centrality
#' @param prepost whether to get presynapses, postsynapses or both
#'
#' @return a \code{data.frame}
#'
#' @examples
#' \donttest{
#' # Choose a bodyid
#' id = c("452677169")
#'
#' # Read in these neurons
#' neuron = neuprint_read_neurons(id)
#'
#' # Re-root
#' neuron.flow = flow_centrality(neuron, polypre = TRUE,
#' mode = "centrifugal",
#' split = "distance")
#'
#' # Let's check that this worked
#' syns = hemibrain_extract_synapses(neuron.flow)
#'
#' \dontrun{
#' # See result
#' library(nat)
#' plot3d_split(neuron.flow)
#' points3d(xyzmatrix(subset(syns,prepost==1)), col = "cyan")
#' points3d(xyzmatrix(subset(syns,prepost==0)), col = "red")
#'
#' }}
#' @export
#' @rdname hemibrain_extract_connections
#' @seealso \code{\link{flow_centrality}}
#' @importFrom nat nlapply is.neuron is.neuronlist
hemibrain_extract_synapses <- function(x,
                                       prepost = c("BOTH","PRE","POST"),
                                       ...){
  prepost = match.arg(prepost)
  if(is.neuronlist(x)){
    syns = nlapply(x,extract_synapses, ...)
    syns = do.call(rbind,syns)
  }else if (is.neuron(x)){
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

#' @export
#' @rdname hemibrain_extract_connections
hemibrain_extract_connections <- function(x,
                                       prepost = c("BOTH","PRE","POST"),
                                       ...){
  prepost = match.arg(prepost)
  if(nat::is.neuronlist(x)){
    syns = nat::nlapply(x,extract_synapses, unitary = TRUE, ...)
    syns = do.call(rbind,syns)
  }else if (nat::is.neuron(x)){
    syns = extract_synapses(x, unitary = TRUE)
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

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`
# hidden
extract_synapses <-function(x, unitary = FALSE){
  syn = x$connectors
  if(!nrow(syn)){
    warning("Neuron ", x$bodyid," has no synapses")
  }
  syn$bodyid = nullToNA(x$bodyid)
  syn$Label = nullToNA(x$d$Label[match(syn$treenode_id,x$d$PointNo)])
  if(unitary){
    syn %>%
      dplyr::filter(Label %in% c(2,3)) %>%
      dplyr::mutate(prepost = dplyr::case_when(
        prepost==0 ~ 1,
        prepost==1  ~ 0
      )) %>% # i.e. switch perspective, presynapses connect to postsynaptic partners
      group_by(bodyid, partner, prepost, Label) %>%
      dplyr::mutate(weight = n()) %>%
      dplyr::distinct(bodyid, partner, prepost, Label, weight) %>%
      dplyr::select(bodyid, partner, prepost, Label, weight) %>%
      as.data.frame() ->
      syn
  }
  syn
}
