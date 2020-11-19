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
#'
#' @inheritParams flow_centrality
#'
#' @param prepost whether to get presynapses, postsynapses or both
#'
#' @return a \code{data.frame}. Depending on which synapse function was called, it can contain the columns:
#'
#' \itemize{
#'
#'   \item{"treenode_id"} { - the position of the node in the SWC-style table found at \code{neuron$d}, where the neuron is the skeleton for \code{bodyid}.}
#'
#'   \item{"connector_id"}{ - the unique ID for a pre/post synapse, as read from neuPrint. If this is not given, you are looking at a connection not a synapse.
#'   In this case \code{count} should be given, which shows the number of synapses in this connection.}
#'
#'   \item{"prepost"}{ - whether the given synapse is a pre-synapse (0, output synapse) or postsynapse (1, input synapse). Alternatively, if a connection is given,
#'   whether this connection is presynaptic to \code{bodyid} (0, \code{bodyid} is target) or postsynaptic (1, \code{bodyid} is source).}
#'
#'   \item{"x"}{ - x coordinate for the root point.}
#'
#'   \item{"y"}{ - y coordinate for the root point.}
#'
#'   \item{"z"}{ - z coordinate for the root point.}
#'
#'   \item{"confidence"}{ - FlyEM's confidence level. The lower the score, the more likely this synapse is an artefact.}
#'
#'   \item{"bodyid"}{ - The neuPrint neuron/body related to the synapse/connection given in each row.}
#'
#'   \item{"partner"}{ - The neuron connecting to \code{bodyid} by the given synapse/connection.}
#'
#'   \item{"pre"}{ - The body ID for the presynaptic (source) neuron.}
#'
#'   \item{"partner"}{ - The body ID for the presynaptic (target) neuron.}
#'
#'   \item{"Label"}{ - The compartment of the \code{bodyid} neuron on which the synapse is placed / which receives/makes the given connection.
#'   See \code{?standardise}.}
#'
#'   \item{"partner.Label"}{ - The compartment of the \code{partner} neuron on which the synapse is placed / which receives/makes the given connection.}
#'
#'   \item{"count"}{ - The number of synapses that make the given connection. Sometimes referred to as 'weight'.}
#'
#'   \item{"norm"}{ - The normalised synapse weight. \code{count} is divided by the total number of inputs that the
#'   target neuron's (\code{post}) compartment (\code{Label}) has. I.e. this normalisation is by total inputs onto a dendrite or axon, not the whole neuron.}
#'
#'   \item{"connection"}{ - The type of compartment-compartment connection specified by this row. The first compartment is the source (pre), the second, the target (post).}
#'
#'}
#'
#' @export
#' @rdname hemibrain_extract_connections
#' @seealso \code{\link{flow_centrality}}
#' @importFrom nat nlapply is.neuron is.neuronlist
#' @examples
#' \dontrun{
#' # Choose bodyids
#' some.pns  = sample(pn.ids, 20)
#'
#' # Read in these neurons
#' neurons = neuprintr::neuprint_read_neurons(some.pns)
#'
#' # Re-root
#' neurons.flow = flow_centrality(neurons,
#' polypre = TRUE,
#' mode = "centrifugal", split = "distance")
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
  x = nat::as.neuronlist(x)
  prepost = match.arg(prepost)
  if("bodyid"%in%colnames(x[,])){
    x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  }else if("skid"%in%colnames(x[,])){
    x = add_field_seq(x,x[,"skid"],field="skid")
  }
  if(is.neuronlist(x)){
    syns = nat::nlapply(x,extract_synapses, unitary = FALSE, ...)
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
  x = add_field_seq(x,names(x),field="bodyid")
  prepost = match.arg(prepost)
  if(nat::is.neuronlist(x)){
    syns = nat::nlapply(x, extract_synapses, unitary = TRUE, ...)
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
  syns$Label = standard_compartments(syns$Label)
  rownames(syns) = 1:nrow(syns)
  syns = syns[order(syns$count, decreasing = TRUE),]
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
  syn[[id]] = nullToNA(as.character(x[[id]]))
  syn[[id]] = gsub(" ","",syn[[id]])
  syn$Label = nullToNA(x$d$Label[match(syn$treenode_id,x$d$PointNo)])
  if(unitary){ # connections, rather than synapses
    syn %>%
      dplyr::mutate(prepost = dplyr::case_when(
        .data$prepost==0 ~ 1,
        .data$prepost==1  ~ 0
      )) %>% # i.e. switch perspective, presynapses connect to postsynaptic partners
      dplyr::group_by(.data[[id]], .data$partner, .data$prepost, .data$Label) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::distinct(.data[[id]], .data$partner, .data$prepost, .data$Label, .data$count) %>%
      dplyr::select(.data[[id]], .data$partner, .data$prepost, .data$Label, .data$count) %>%
      as.data.frame(stringsAsFactors = FALSE) ->
      syn
  }
  syn$Label = standard_compartments(syn$Label)
  syn
}

#' @export
#' @rdname hemibrain_extract_connections
hemibrain_extract_compartment_edgelist <- function(x, ...){
  x = add_field_seq(x,names(x),field="bodyid")
  if(nat::is.neuronlist(x)){
    syns.list = nat::nlapply(x, extract_synapses, unitary = FALSE, ...)
  }else{
    stop("x must be a neuronlist object")
  }
  names(syns.list) = NULL
  lookup = nat::nlapply(syns.list, extract_lookup, ...)
  lookup = unlist(lookup)
  sel.cols = c("total.outputs", "total.inputs", "axon.outputs",
               "dend.outputs", "axon.inputs", "dend.inputs")
  if(sum(sel.cols%in%colnames(x[,]))!=length(sel.cols)){
    mets = hemibrain_compartment_metrics(x, OmitFailures = TRUE, delta = 5, resample = NULL, locality = FALSE)
    mets[is.na(mets)] = 0
  }else{
    mets = x[,]
  }
  comp.meta = mets[,c("total.outputs", "total.inputs", "axon.outputs",
                   "dend.outputs", "axon.inputs", "dend.inputs")]
  elists = nat::nlapply(syns.list, extract_elist, lookup = lookup, meta = comp.meta, ...)
  elist = do.call(rbind, elists)
  rownames(elist) = 1:nrow(elist)
  elist = elist[order(elist$norm, decreasing = TRUE),]
  elist = elist[order(elist$count, decreasing = TRUE),]
  elist
}

# hidden, for one pre neuron
extract_elist <- function(syns, lookup, meta = NULL){
  bodyid = unique(nullToNA(as.character(syns$bodyid)))
  if(!is.null(meta)){
    if(!bodyid%in%rownames(meta)){
      stop(bodyid, " not in rownames of meta")
    }
    d.post = nullToNA(as.numeric(meta[bodyid,"dend.inputs"]))
    a.post = nullToNA(as.numeric(meta[bodyid,"axon.inputs"]))
  }else{
    d.post = nrow(subset(syns, syns$prepost==1 & syns$Label==3))
    a.post = nrow(subset(syns, syns$prepost==1 & syns$Label==2))
    warning("Argument 'meta' not given. norm column represents normalisation by the total number of pre-post connections
    in the final data.frame, not necessarily the total number of postsynapses for
    the 'post' neuron.")
  }
  syns %>%
    # Re-name for clarity
    dplyr::filter(.data$prepost==1) %>%
    dplyr::rename(post = .data$bodyid) %>%
    dplyr::rename(pre = .data$partner) %>%
    dplyr::rename(post.Label = .data$Label) %>%
    # Compartment labels
    dplyr::mutate(pre.Label = lookup[as.character(.data$connector_id)]) %>%
    dplyr::mutate(pre.Label = ifelse(is.na(.data$pre.Label),"error",.data$pre.Label)) %>%
    # Synapse counts
    dplyr::group_by(.data$post, .data$pre, .data$post.Label, .data$pre.Label) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    # Normalised synapses, by compartment
    dplyr::ungroup() %>%
    dplyr::group_by(.data$post,.data$post.Label) %>%
    dplyr::mutate(norm = .data$count/ifelse(.data$post.Label%in%c(3,"dendrite","dendrites","dend"),d.post,a.post)) %>%
    # Clean up
    dplyr::distinct(.data$post, .data$pre,.data$post.Label, .data$pre.Label, .data$count, .data$norm) %>%
    dplyr::select(.data$post, .data$pre, .data$post.Label, .data$pre.Label, .data$count, .data$norm) %>%
    dplyr::filter(!is.na(.data$pre.Label) & .data$count > 0) %>%
    as.data.frame(stringsAsFactors = FALSE) ->
    elist
  rownames(elist) = 1:nrow(elist)
  elist$post.Label = standard_compartments(elist$post.Label)
  elist$pre.Label = standard_compartments(elist$pre.Label)
  elist$connection = paste(elist$pre.Label,elist$post.Label,sep="-")
  subset(elist, elist$pre.Label!="error")
}

# hidden
extract_lookup <- function(syns){
  syns %>%
    dplyr::filter(.data$prepost==0) %>%
    dplyr::distinct(.data$bodyid, .data$partner, .data$connector_id, .data$Label) %>%
    as.data.frame(stringsAsFactors = FALSE) ->
    conn.lookup
  lookup = conn.lookup$Label
  names(lookup) = conn.lookup$connector_id
  lookup
}
