#' Flexible grouping of rows/columns of a raw adjacency matrix
#'
#' Use this to create a summary adjacency matrix in which input and/or output
#' neurons are coalesced into a smaller number of rows/columns. This can be used
#' as a preprocessing step to present a simple representation of the
#' connectivity between a larger number of cells in a heatmap.
#'
#' @param inputids Either the bodyids of the input neurons OR an adjacency
#'   matrix
#' @param outputids bodyids of the output neurons. Not required if
#'   \code{inputids} is an adjacency matrix.
#' @param ingroup,outgroup grouping variables used to coalesce related neurons
#' @param threshold Remove (groups of) neurons that make fewer connections than
#'   this.
#' @param scale Whether to scale the rows or columns so that they sum to 1
#' @param ... Additional arguments passed to
#'   \code{\link{neuprint_get_adjacency_matrix}}
#'
#' @return named numeric matrix
#' @export
#'
#' @examples
#' \donttest{
#' pnkc=grouped_adjacency(class2ids("PN"), 'KC')
#' heatmap(pnkc)
#' }
#' @importFrom neuprintr neuprint_get_adjacency_matrix neuprint_get_meta
grouped_adjacency <- function(inputids=NULL, outputids=NULL,
                              ingroup=c("type", "name", "cellBodyFiber"),
                              outgroup=c("type", "name", "cellBodyFiber"),
                              threshold=0,
                              scale=c("none", "col", "row"), ...) {

  if(is.matrix(inputids)){
    am <- inputids
    inputids=rownames(am)
    outputids=colnames(am)
  } else am <- NULL

  if(is.null(inputids)) {
    inputids=names(ingroup)
    if(is.null(inputids))
      stop("You must supply valid input ids!")
  }
  inputids <- neuprint_ids(inputids)
  if(is.function(ingroup))
    ingroup = ingroup(inputids)
  else if(is.character(ingroup)) {
    ingroup=match.arg(ingroup)
    inm = neuprint_get_meta(inputids)
    # in case some neuron doesn't have metadata
    ingroup <- as.character(inm[[ingroup]][match(inputids, inm$bodyid)])
  } else if((is.numeric(ingroup)|| is.factor(ingroup)) && !is.null(names(ingroup))) {
    # input looks numeric and has names
    stopifnot(all(names(ingroup) %in% inputids))
    inputids=intersect(names(ingroup), inputids)
  }

  if(is.null(outputids)) {
    outputids=names(outgroup)
    if(is.null(outputids))
      stop("You must supply valid input ids!")
  }
  outputids <- neuprint_ids(outputids)
  outm <- NULL
  if(is.function(outgroup))
    outgroup = outgroup(outputids)
  else if(is.character(outgroup)) {
    outgroup=match.arg(outgroup)
    outm = neuprint_get_meta(outputids)
    # in case some neuron doesn't have metadata
    outgroup <- as.character(outm[[outgroup]][match(outputids, outm$bodyid)])
  } else if((is.numeric(outgroup)|| is.factor(outgroup)) && !is.null(names(outgroup))) {
    # input looks numeric and has names
    stopifnot(all(names(outgroup) %in% outputids))
    outputids=intersect(names(outgroup), outputids)
  }

  if(is.null(am))
    am <- neuprint_get_adjacency_matrix(inputids = inputids, outputids = outputids, ...)
  else {
    # in case the grouping variable implied a reordering
    am <- am[inputids, outputids, drop=FALSE]
  }

  if(!is.null(ingroup)) {
    l=split(as.data.frame(am), ingroup)
    # nb transpose back to original orientation
    am=t(sapply(l, colSums))
  }
  if(!is.null(outgroup)) {
    l=split(as.data.frame(t(am)), outgroup)
    # don't transpose so back to inputs as rows
    am=sapply(l, colSums)
  }

  scale=match.arg(scale)
  if(scale=="row") {
    am <- t(scale(t(am), center = FALSE, scale=rowSums(am)))
  } else if(scale=='col') {
    am <- scale(am, center = FALSE, scale=colSums(am))
  }
  if(threshold>0) am=am[, !colSums(am)<threshold, drop=F]
  if(threshold>0) am=am[!rowSums(am)<threshold, , drop=F]
  as.matrix(am)
}

