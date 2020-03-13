############################################################################
################################ Synapses ##################################
############################################################################

#' Extract synapse location and labels, and edgelists from a neuron/neuronlist
#'
#' @description Extract a single data frame describing synapse/connection types,
#'   partners, locations and position on a neuron's axon/dendrite. You can
#'   either get all synapses returned or all unitary connections to a neuron's
#'   partners returned. Broken down by axon/dendrite (\code{Label}), and
#'   pre/postsynapses or pre/postsynaptic partners. Note that \code{hemibrain_extract_compartment_edgelist} will
#'   only return connections between neurons given in the argument \code{x}.
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
    syns = nlapply(x,extract_synapses, unitary = FALSE, ...)
    syns = do.call(rbind,syns)
  }else if (is.neuron(x)){
    syns = extract_synapses(x, unitary = FALSE)
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
  x = add_field_seq(x,x[,"bodyid"],field="bodyid")
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
#' @importFrom dplyr filter mutate group_by distinct select n case_when
#' @importFrom rlang .data
extract_synapses <-function(x, unitary = FALSE){
  syn = x$connectors
  if(!nrow(syn)){
    warning("Neuron ", x$bodyid," has no synapses")
  }
  syn$bodyid = nullToNA(x$bodyid)
  syn$Label = nullToNA(x$d$Label[match(syn$treenode_id,x$d$PointNo)])
  if(unitary){ # connections, rather than synapses
    syn %>%
      dplyr::mutate(prepost = dplyr::case_when(
        .data$prepost==0 ~ 1,
        .data$prepost==1  ~ 0
      )) %>% # i.e. switch perspective, presynapses connect to postsynaptic partners
      dplyr::group_by(.data$bodyid, .data$partner, .data$prepost, .data$Label) %>%
      dplyr::mutate(weight = dplyr::n()) %>%
      dplyr::distinct(.data$bodyid, .data$partner, .data$prepost, .data$Label, .data$weight) %>%
      dplyr::select(.data$bodyid, .data$partner, .data$prepost, .data$Label, .data$weight) %>%
      as.data.frame() ->
      syn
  }
  syn
}

#' @export
#' @rdname hemibrain_extract_connections
hemibrain_extract_compartment_edgelist <- function(x, ...){
  x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  if(nat::is.neuronlist(x)){
    syns = nat::nlapply(x, extract_synapses, unitary = FALSE, ...)
    syns = do.call(rbind,syns)
  }else if (nat::is.neuron(x)){
    syns = extract_synapses(x, unitary = FALSE)
  }else{
    stop("x must be a neuron or neuronlist object")
  }
  syns %>%
    dplyr::filter(prepost==1) %>%
    dplyr::distinct(.data$bodyid, .data$partner, .data$connector_id, .data$Label) %>%
    as.data.frame() ->
    conn.lookup
  lookup = conn.lookup$Label
  names(lookup) = conn.lookup$connector_id
  syns %>%
    dplyr::filter(prepost==0) %>%
    dplyr::mutate(partner.Label = lookup[as.character(connector_id)]) %>%
    dplyr::group_by(.data$bodyid, .data$partner, .data$Label) %>%
    dplyr::mutate(weight = dplyr::n()) %>%
    dplyr::distinct(.data$bodyid, .data$partner,.data$Label, .data$partner.Label, .data$weight) %>%
    dplyr::select(.data$bodyid, .data$partner, .data$Label, .data$partner.Label, .data$weight) %>%
    dplyr::filter(!is.na(partner.Label)) %>%
    as.data.frame() ->
    elist
  rownames(elist) = 1:nrow(elist)
  elist$Label = standard_compartments(elist$Label)
  elist$partner.Label = standard_compartments(elist$partner.Label)
  elist
}

