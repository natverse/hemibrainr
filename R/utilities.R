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


# hidden
maxout <- function(x, max){
  x[x>max] = max
  x[is.infinite(x)|is.null(x)|is.na(x)] = 0
  x
}

# hidden
lengthnorm <- function(x){
  sum(x, na.rm = TRUE)/length(x)
}

# hidden
prune_synapseless_branches <- function(x){
  s = x$SubTrees
  prune = c()
  if(is.null(x$connectors$treenode_id)){
    stop("No connectors in neuron,")
  }
  for(t in 1:x$nTrees){
    ss = unlist(s[[t]])
    with.syns = sum(ss%in%x$connectors$treenode_id)>0
    if(!with.syns){
      prune = c(prune, ss)
    }
  }
  if(length(prune)>0){
    y = nat::prune_vertices(x, verticestoprune = prune, invert = FALSE)
    y$connectors = x$connectors[x$connectors$treenode_id %in% y$d$PointNo, ]
    relevant.points = subset(x$d, PointNo %in% y$d$PointNo)
    y$d = relevant.points[match(y$d$PointNo, relevant.points$PointNo), ]
  }else{
    y=x
  }
  y = hemibrain_neuron_class(y)
  y
}

# hidden
hemibrain_neuron_class <- function (x){
  class(x) = unique(c(class(x),"neuprintneuron","catmaidneuron","neuron","list"))
  x
}
