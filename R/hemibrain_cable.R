###################################################################################
################################ Get neuron cable #################################
###################################################################################

#' Retrieve the primary neurite (cell body fibre) for a neuron
#'
#' @description Extract the cell body fibre for an insect neuron,
#' assuming that this neuron has a soma and that the soma is the neuron's root.
#'
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param neuron logical, whether to return a neuron/neuronlist object (TRUE) or
#' the IDs for nodes in the primary neurite
#' @param invert logical, whether to return the primary neurite (TRUE) or the rest
#' of the neuron with its primary neurite removed (FALSE).
#' (note: not x$d$PointNo, but rownames(x$d)).
#' @param ... methods sent to \code{nat::nlapply}
#'
#' @return the neuron or neuron list object inputted, with centripetal flow
#'   centrality information added to neuron$d and a segregation index score.
#'   The neuron$d$Label now gives the compartment, where axon is Label = 2,
#'   dendrite Label = 3, primary dendrite Label = 9 and primary neurite Label = 7.
#'   Soma is Label = 1.
#'
#' @examples
#' \donttest{
#' # Get neuron
#' neuron = neuprintr::neuprint_read_neurons("451987038")
#'
#' # Extract primary neurite
#' pnt = primary_neurite(neuron)
#'
#' \dontrun{
#' # Plot the primary neurite
#' nat::nopen3d()
#' plot3d(pnt, lwd = 2, col = "#C70E7B", soma = 1000)
#' plot3d(neuron, lwd = 2, col = "grey30")
#' }}
#' @export
#' @seealso \code{\link{flow_centrality}}
primary_neurite <-function(x, neuron = TRUE, invert = FALSE, ...) UseMethod("primary_neurite")

#' @export
primary_neurite.neuron <- function(x, neuron = TRUE, invert = FALSE, ...){
  simp = nat::simplify_neuron(x, n=1)
  pnt = sapply(simp$SegList,function(s) sum(simp$StartPoint%in%s)>0)
  pnt = simp$SegList[[which(pnt&!duplicated(pnt))]]
  if(neuron){
    pnt = nat::prune_vertices(simp, verticestoprune = pnt, invert = !invert)
  }else{
    pnt = (1:nrow(x$d))[match(simp$d$PointNo[pnt],x$d$PointNo)]
    if(invert){
      pnt = setdiff(1:nrow(x$d),pnt)
    }
  }
  pnt
}

#' @export
primary_neurite.neuprintneuron <- primary_neurite.neuron

#' @export
primary_neurite.catmaidneuron <- primary_neurite.neuron

#' @export
primary_neurite.neuronlist <- function(x, neuron = TRUE, invert = FALSE, ...){
  nat::nlapply(x, primary_neurite.neuron, neuron = neuron, invert = invert,  ...)
}


#' Extract axonic/dendritic points/cable from a neuron/neuronlist
#'
#' @description Extract axonic/dendritic/mixed/primary dendrite points/endpoints/cable from a neuron/neuronlist object
#'
#' @param x a neuron/neuronlist object that has its axons/dendrites labelled in swc format in its neuron$d dataframes
#' @param mixed whether or not to include points assigned as uncertain or mixed polarity cable
#' @param ... additional arguments passed to methods
#'
#' @return a matrix of 3D points
#' @export
#' @rdname extract_cable
axonic_points<-function(x, ...) UseMethod("axonic_points")
#' @export
#' @rdname extract_cable
dendritic_points<-function(x, ...) UseMethod("dendritic_points")
#' @export
#' @rdname extract_cable
mixed_points<-function(x, ...) UseMethod("mixed_points")
#' @export
#' @rdname extract_cable
primary_dendrite_points<-function(x, ...) UseMethod("primary_dendrite_points")

#' @export
axonic_points.neuron <- function(x, ...){
  points=x$d
  nat::xyzmatrix(points[points$Label%in%c(-2,2),])
}

#' @export
dendritic_points.neuron <- function(x, ...){
  points=x$d
  nat::xyzmatrix(points[points$Label%in%c(-3,3),])
}

#' @export
mixed_points.neuron <- function(x, ...){ # Mised also means that I do not know
  points=x$d
  nat::xyzmatrix(points[points$Label%in%c(8),])
}

#' @export
primary_dendrite_points.neuron <- function(x, ...){ # Mised also means that I do not know
  points=x$d
  nat::xyzmatrix(points[points$Label%in%c(4),])
}

#' @export
dendritic_points.neuronlist <- function(x, ...){
  do.call(rbind,nat::nlapply(x,dendritic_points.neuron, ...))
}

#' @export
axonic_points.neuronlist <- function(x, ...){
  do.call(rbind,nat::nlapply(x,axonic_points.neuron, ...))
}

#' @export
mixed_points.neuronlist <- function(x, ...){
  do.call(rbind,nat::nlapply(x, mixed_points.neuron, ...))
}

#' @export
primary_dendrite_points.neuronlist <- function(x, ...){
  do.call(rbind,nat::nlapply(x, primary_dendrite_points.neuron, ...))
}

#' @export
#' @rdname extract_cable
axonal_endings <- function(x){
  points=x$d[nat::endpoints(x)[which(nat::endpoints(x)!=nat::rootpoints(x))],]
  nat::xyzmatrix(points[points$Label%in%c(-2,2),])
}

#' @export
#' @rdname extract_cable
dendritic_endings <- function(x){
  points=x$d[nat::endpoints(x)[which(nat::endpoints(x)!=nat::rootpoints(x))],]
  nat::xyzmatrix(points[points$Label%in%c(-3,3),])
}

#' @export
#' @rdname extract_cable
axonic_endings <- function(x){
  points=x$d[nat::endpoints(x)[which(nat::endpoints(x)!=nat::rootpoints(x))],]
  nat::xyzmatrix(points[points$Label%in%c(-2,2),])
}

#' @export
#' @rdname extract_cable
primary_dendrite_endings <- function(x){
  if(nat::is.neuron(x)){
    x = primary_dendrite_cable.neuron(x)
    points=x$d[nat::endpoints(x),]
  }else{
    nat::nlapply(x,function(x) primary_dendrite_cable.neuron(x)$d[nat::endpoints(primary_dendrite_cable.neuron(x)),])
  }
}

#' @export
#' @rdname extract_cable
axonic_cable<-function(x, mixed = FALSE, ...) UseMethod("axonic_cable")

#' @export
#' @rdname extract_cable
dendritic_cable<-function(x, mixed = FALSE, ...) UseMethod("dendritic_cable")

#' @export
#' @rdname extract_cable
arbour_cable<-function(x, mixed = FALSE, ...) UseMethod("arbour_cable")

#' @export
#' @rdname extract_cable
unsure_cable<-function(x, mixed = FALSE, ...) UseMethod("unsure_cable")

#' @export
#' @rdname extract_cable
primary_dendrite_cable<-function(x, ...) UseMethod("primary_dendrite_cable")

#' @export
axonic_cable.neuron <- function(x, mixed=FALSE, ...){
  points=x$d
  if (mixed==TRUE){
    chosen = c(-2,2,8)
  }else{
    chosen = c(-2,2)
  }
  v = subset(rownames(x$d), x$d$Label %in% chosen)
  if("neuprintneuron"%in%class(x)){
    neuron = prune_vertices.neuprintneuron(x=x,verticestoprune=v,invert=TRUE)
  }else{
    neuron = nat::prune_vertices(x,verticestoprune=v,invert=TRUE)
  }
  neuron$d$Label=2
  neuron
}

#' @export
axonic_cable.neuprintneuron <- axonic_cable.neuron

#' @export
dendritic_cable.neuron <- function(x, mixed = FALSE, ...){
  points=x$d
  if (mixed==T){
    chosen = c(-3,3,8)
  } else{
    chosen = c(-3,3)
  }
  v = subset(rownames(x$d), x$d$Label %in% chosen)
  if("neuprintneuron"%in%class(x)){
    neuron = prune_vertices.neuprintneuron(x,verticestoprune=v,invert=TRUE)
  }else{
    neuron = nat::prune_vertices(x,verticestoprune=v,invert=TRUE)
  }
  neuron$d$Label=3
  neuron
}

#' @export
dendritic_cable.neuprintneuron <- dendritic_cable.neuron

#' @export
arbour_cable.neuron <- function(x, mixed = FALSE, ...){
  points=x$d
  if (mixed==T){
    chosen = c(-3,3,2,-2,8)
  }else{
    chosen = c(-3,3,2,-2)
  }
  v = subset(rownames(x$d), x$d$Label %in% chosen)
  if("neuprintneuron"%in%class(x)){
    neuron = prune_vertices.neuprintneuron(x,verticestoprune=v,invert=TRUE)
    class(neuron) = c("neuprintneuron","neuron")
  }else{
    neuron = nat::prune_vertices(x,verticestoprune=v,invert=TRUE)
  }
  neuron
}


#' @export
arbour_cable.neuprintneuron <- arbour_cable.neuron

#' @export
unsure_cable.neuron <- function(x, mixed=FALSE, ...){
  points=x$d
  chosen = c(-8,8:100)
  v = subset(rownames(x$d), x$d$Label %in% chosen)
  if("neuprintneuron"%in%class(x)){
    neuron = prune_vertices.neuprintneuron(x,verticestoprune=v,invert=TRUE)
    class(neuron) = c("neuprintneuron","neuron")
  }else{
    neuron = nat::prune_vertices(x,verticestoprune=v,invert=TRUE)
  }
  neuron$d$Label=8
  neuron
}

#' @export
unsure_cable.neuprintneuron <- unsure_cable.neuron

#' @export
primary_dendrite_cable.neuron <- function(x, ...){
  points=x$d
  v = subset(rownames(x$d), x$d$Label %in% 4)
  if("neuprintneuron"%in%class(x)){
    neuron = prune_vertices.neuprintneuron(x,verticestoprune=v,invert=TRUE)
  }else{
    neuron = nat::prune_vertices(x,verticestoprune=v,invert=TRUE)
  }
  neuron$d$Label=3
  neuron
}

#' @export
primary_dendrite_cable.neuprintneuron <- primary_dendrite_cable.neuron

#' @export
axonic_cable.neuronlist <- function(x, mixed=FALSE, ...){
  nat::nlapply(x,axonic_cable.neuron,mixed=mixed,OmitFailures = T, ...)
}

#' @export
dendritic_cable.neuronlist <- function(x, mixed=FALSE, ...){
  nat::nlapply(x,dendritic_cable.neuron,mixed=mixed,OmitFailures = T, ...)
}

#' @export
arbour_cable.neuronlist <- function(x, mixed=FALSE, ...){
  nat::nlapply(x,arbour_cable.neuron,mixed=mixed,OmitFailures = T, ...)
}

#' @export
unsure_cable.neuronlist <- function(x, mixed = FALSE, ...){
  nat::nlapply(x,unsure_cable.neuron,OmitFailures = T, ...)
}

#' @export
primary_dendrite_cable.neuronlist <- function(x, ...){
  nat::nlapply(x,primary_dendrite_cable.neuron,OmitFailures = T, ...)
}

#' Prune vertices from a neuprint neuron (taking care of synapse etc information)
#'
#' @export
#' @importFrom nat prune_vertices
#' @inheritParams nat::prune_vertices
#' @seealso \code{nat::\link[nat]{prune_vertices}}
prune_vertices.neuprintneuron <- function (x, verticestoprune, invert = FALSE, ...){
  pruned = nat::prune_vertices(x, verticestoprune, invert = invert,
                               ...)
  pruned$connectors = x$connectors[x$connectors$treenode_id %in%
                                     pruned$d$PointNo, ]
  relevant.points = subset(x$d, x$d$PointNo %in% pruned$d$PointNo)
  y = pruned
  y$d = relevant.points[match(pruned$d$PointNo, relevant.points$PointNo),
                        ]
  y$d$Parent = pruned$d$Parent
  y$tags = lapply(x$tags, function(t) t[t %in% pruned$d$PointNo])
  y$url = x$url
  y$headers = x$headers
  y$AD.segregation.index = x$AD.segregation.index
  class(y) = c("neuprintneuron", "neuron")
  y
}
