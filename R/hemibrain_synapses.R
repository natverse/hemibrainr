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
#' @export
#' @rdname hemibrain_extract_connections
#' @seealso \code{\link{flow_centrality}}
#' @importFrom nat nlapply is.neuron is.neuronlist
#' @examples
#' \dontrun{
#' # Choose bodyids
#' al.local.neurons =
#'  c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#'  "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#'  "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#' "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#'  "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#'  "1794037618", "5813018729", "2013333009")
#'
#' # Read in these neurons
#' neurons = neuprintr::neuprint_read_neurons(al.local.neurons)
#'
#' # Re-root
#' neurons.flow = flow_centrality(neurons, polypre = TRUE,
#'   mode = "centrifugal", split = "distance")
#'
#' # Let's check that this worked
#' syns = hemibrain_extract_synapses(neurons.flow)
#'
#' # Get the edgelist by compartment
#' elist = hemibrain_extract_compartment_edgelist(neurons.flow)
#'
#' # See result
#' library(nat)
#' plot3d_split(neuron.flow)
#' points3d(xyzmatrix(subset(syns,prepost==1)), col = "cyan")
#' points3d(xyzmatrix(subset(syns,prepost==0)), col = "red")
#' }
hemibrain_extract_synapses <- function(x,
                                       prepost = c("BOTH","PRE","POST"),
                                       ...){
  prepost = match.arg(prepost)
  if("bodyid"%in%colnames(x[,])){
    x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  }else if("skid"%in%colnames(x[,])){
    x = add_field_seq(x,x[,"skid"],field="skid")
  }
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
  if(!is.null(x$bodyid)){
    id = "bodyid"
  }else if(!is.null(x$skid)){
    id = "skid"
  }else{
    id = "id"
  }
  syn[[id]] = nullToNA(x[[id]])
  syn[[id]] = gsub(" ","",syn[[id]])
  syn$Label = nullToNA(x$d$Label[match(syn$treenode_id,x$d$PointNo)])
  if(unitary){ # connections, rather than synapses
    syn %>%
      dplyr::mutate(prepost = dplyr::case_when(
        .data$prepost==0 ~ 1,
        .data$prepost==1  ~ 0
      )) %>% # i.e. switch perspective, presynapses connect to postsynaptic partners
      dplyr::group_by(.data[[id]], .data$partner, .data$prepost, .data$Label) %>%
      dplyr::mutate(weight = dplyr::n()) %>%
      dplyr::distinct(.data[[id]], .data$partner, .data$prepost, .data$Label, .data$weight) %>%
      dplyr::select(.data[[id]], .data$partner, .data$prepost, .data$Label, .data$weight) %>%
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
    syns.list = nat::nlapply(x, extract_synapses, unitary = FALSE, ...)
  }else{
    stop("x must be a neuronlist object")
  }
  names(syns.list) = NULL
  lookup = nat::nlapply(syns.list, extract_lookup, ...)
  lookup = unlist(lookup)
  elists = nat::nlapply(syns.list, extract_elist, lookup = lookup, ...)
  elist = do.call(rbind, elists)
  rownames(elist) = 1:nrow(elist)
  elist
}

# hidden
extract_elist <- function(syns, lookup){
  syns %>%
    dplyr::filter(.data$prepost==0) %>%
    dplyr::mutate(partner.Label = lookup[as.character(.data$connector_id)]) %>%
    dplyr::group_by(.data$bodyid, .data$partner, .data$Label) %>%
    dplyr::mutate(weight = dplyr::n()) %>%
    dplyr::distinct(.data$bodyid, .data$partner,.data$Label, .data$partner.Label, .data$weight) %>%
    dplyr::select(.data$bodyid, .data$partner, .data$Label, .data$partner.Label, .data$weight) %>%
    dplyr::filter(!is.na(.data$partner.Label) & .data$weight > 0) %>%
    as.data.frame() ->
    elist
  rownames(elist) = 1:nrow(elist)
  elist$Label = standard_compartments(elist$Label)
  elist$partner.Label = standard_compartments(elist$partner.Label)
  elist
}

# hidden
extract_lookup <- function(syns){
  syns %>%
    dplyr::filter(.data$prepost==1) %>%
    dplyr::distinct(.data$bodyid, .data$partner, .data$connector_id, .data$Label) %>%
    as.data.frame() ->
    conn.lookup
  lookup = conn.lookup$Label
  names(lookup) = conn.lookup$connector_id
  lookup
}
