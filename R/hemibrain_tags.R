###############################################################################
################################ Tag neurons ##################################
###############################################################################

#' Set and see tags on a neuron object inc. within a neuronlist
#'
#' @description Set the x$tags field within a neuron object.
#'
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object.
#' @param y a neuron object, which is to receive tags from x, located at \code{x$tags}.
#' @param i the index of the arguments provides through ..., to be used as the tag, if these arguments
#' are a vector of greater than length 1. Principally for use within method \code{hemibrainr:::list}.
#' @param ... named arguments with values, to be added to x$tags. E.e. using soma = TRUE,
#' will add \code{x$tags$soma}, where \code{x$tags$soma} is \code{TRUE}. Tags are expected to be of length 1.
#' If x is a \code{neuronlist} then tags may be supplied of \code{length(x)}.
#'
#' @examples
#' \donttest{
#' # Choose neurons
#' ## In this case some antennal lobe local neurons
#' al.local.neurons = c("1702323386", "2068966051")
#'
#' # Get neurons
#' neurons = hemibrain_read_neurons(al.local.neurons)
#'
#' # Tag
#' neurons = hemibrain_settags(neurons, antennal.lobe = TRUE)
#'
#' # See tags
#' hemibrain_seetags(x)
#' }
#' @rdname hemibrain_tags
#' @export
#' @seealso \code{\link{hemibrain_read_neurons}}
hemibrain_settags <-function(x, i = 1, ...) UseMethod("hemibrain_settags")
#' @rdname hemibrain_tags
#' @export
hemibrain_settags.neuron <- function(x, i = 1, ...){
  m = as.list(match.call())
  tags = as.character(setdiff(names(m), c("hemibrain_settags","x","i","")))
  m = lapply(tags, function(t) do.call(paste0, list(m[[t]], sep="", collapse=" ... ")))
  m =  lapply(m, strsplit, split = " ... ")
  names(m) = tags
  for(tag in tags){
    t = unname(unlist(nullToNA(m[[tag]])))[i]
    t = ifelse(is.na(t),"none",t)
    x$tags[[tag]] = t
  }
  x
}
#' @rdname hemibrain_tags
#' @export
hemibrain_settags.neuronlist<-function(x, ...){
  z = lapply(seq_along(x), function(y) hemibrain_settags.neuron(x[[y]], i = y, ...))
  z = nat::as.neuronlist(z)
  names(z) = names(x)
  z[,] = x[,]
  z
}
#' @rdname hemibrain_tags
#' @export
hemibrain_seetags <-function(x, ...) UseMethod("hemibrain_seetags")
#' @rdname hemibrain_tags
#' @export
hemibrain_seetags.neuron <- function(x,...){
  bodyid = x$bodyid
  sop = data.frame(stringsAsFactors = FALSE)
  for(tag in names(x$tags)){
    col = nullToNA(x$tags[[tag]])
    sop = rbind(sop,col)
  }
  sop = t(sop)
  colnames(sop) = names(x$tags)
  sop = cbind(sop, data.frame(bodyid = bodyid, stringsAsFactors = FALSE))
  rownames(sop) = NULL
  as.data.frame(sop, stringsAsFactors = FALSE)
}
#' @rdname hemibrain_tags
#' @export
hemibrain_seetags.neuronlist<-function(x,...){
  data = lapply(x,hemibrain_seetags.neuron,...)
  data = do.call(plyr::rbind.fill, data)
  data[is.na(data)] = "none"
  data
}
#' @rdname hemibrain_tags
#' @export
hemibrain_carryover_tags <- function(x, y){
  if(!nat::is.neuron(x)){
    warning("x must be a neuron object")
    return(x)
  }
  tags = names(x$tags)
  if(!is.issue(tags)){
    for(t in tags){
      y$tags[t] <- x$tags[[t]]
    }
  }
  y$bodyid = x$bodyid
  y$soma = x$soma
  y
}
#' @rdname hemibrain_tags
#' @export
metadata_add_tags <- function(x){
  if(!nat::is.neuronlist(x)){
    warning("x must be a neuronlist object")
    return(x)
  }
  x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  tags = hemibrain_seetags(x)
  tags$somapos = tags$soma
  tags$soma = NULL
  df = (x[,c( "bodyid",setdiff( colnames(x[,]), colnames(tags)) )])
  df = merge(df, tags, all.x = TRUE, all.y = FALSE)
  rownames(df) = df$bodyid
  x[,] = df[names(x),]
  x
}

#' Add a field to neuron objects (inc. in a neuronlist)
#'
#' @description Add an item to the the list object that comprises a neuron
#'
#' @inheritParams flow_centrality
#' @param entry item to add.
#' @param entries a vector/list of entries, the same length as \code{x}.
#' Each neuron in \code{x} will have the entry in entries as the same index, added to it
#' as \code{field}.
#' @param field name of new field, or field to overwrite.
#'
#' @return a \code{neuron} or \code{neuronlist}
#' @seealso \code{\link{add_Label}}
#' @export
add_field <-function(x, entry, field = "bodyid", ...) UseMethod("add_field")

#' @export
add_field.neuron <- function(x, entry, field = "bodyid", ...){
  x[[field]] = entry
  x
}
#' @export
add_field.neuronlist <- function(x, entry, field = "bodyid", ...){
  nat::nlapply(x, add_Label.neuron, entry, field = "bodyid", ...)
}
#' @export
#' @rdname add_field
add_field_seq <- function(x, entries, field = "bodyid", ...){
  x = nat::as.neuronlist(x)
  if(length(entries)!=length(x)){
    stop("The length of the entries to add must be the same as the length of the neuronlist, x")
  }
  nl = nat::neuronlist()
  for(i in 1:length(x)){
    y = x[[i]]
    entry = entries[i]
    y = add_field(y, entry = entry, field = field, ...)
    nl = c(nl, nat::as.neuronlist(y))
  }
  names(nl) = names(x)
  nl[,] = x[,]
  nl
}
