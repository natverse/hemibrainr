#############################################################################
################################ Utilities ##################################
#############################################################################

# hidden
nullToNA <- function(x) {
  if(is.list(x)){
    x[sapply(x, is.null)] <- NA
  }else{
    x = sapply(x, function(y) ifelse(is.null(y)|!length(y), NA, y))
    if(!length(x)){
      x = NA
    }
  }
  x
}

# hidden
break_into_subtrees <- function(x){
  nlist = nat::neuronlist()
  for(i in 1:x$nTrees){
    segs = x$SubTrees[[i]]
    d = x$d[unique(unlist(segs)),]
    rownames(d) = 1:nrow(d)
    n = nat::as.neuron(d)
    nlist  = c(nlist, nat::as.neuronlist(n))
  }
  nlist
}

# hidden
is.hxsurf <- function(x){
  "hxsurf"%in%class(x)
}

#' @importFrom nat progress_natprogress
#' @export
nat::progress_natprogress

# Get the Nth highest number
maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}
