# function to do with neuron morphology

#' Calculate a connectivity similarity score between two connectivity profiles
#'
#' @description Calculate a similarity score between connectivity matrices that penalises small differences between low and does not heavily penalise large differences between high weights. Algorithm from Jarrell et al. 2012.
#' @param x a \code{neuron} or \code{neuronlist} object.
#' @param Label defaults to 1, so distances to soma are chiosen. If 7, then distance from primary neurite's branchpoint.
#' @param graph.distance logical, whether to use geodesic distance in terms of node number (\code{TRUE}) or real unit distance.
#' @param ... methods passed to \code{nat::nlapply}
#' @name geodesic_distance
geodesic_distance <- function(x, graph.distance=TRUE, Label = 1, ...) UseMethod("geodesic_distance")

#' @export
#' @rdname geodesic_distance
geodesic_distance.neuron <- function(x, graph.distance=TRUE, Label = 1, ...){
  if(Label == 1){
    bp = as.numeric(nat::rootpoints(x))
  }else{
    pn = subset(x$d, Label == Label)
    bp = nat::endpoints(pn)[!nat::endpoints(pn) %in% nat::rootpoints(pn)]
    bp = as.numeric(rownames(subset(x$d, PointNo == bp)))
  }
    n = nat::as.ngraph(x)
    path = suppressWarnings(igraph::shortest_paths(n, from = bp,
                                                   mode = "out")$vpath)
    x$d$geodesic.distance = sapply(path, length)
    if (!graph.distance) {
      conns = as.numeric(rownames((subset(x$d, PointNo %in%
                                            x$connectors$treenode_id))))
      paths = suppressWarnings(igraph::shortest_paths(n, from = bp,
                                                      to = conns, mode = "out")$vpath)
      real.lengths = c()
      for (p in paths) {
        if (length(p) > 2) {
          rl = summary(nat::prune_vertices(x, verticestoprune = p,
                                           invert = TRUE))$cable.length
          real.lengths = c(real.lengths, rl)
        }
        else {
          real.lengths = c(real.lengths, 0)
        }
      }
      x$d$geodesic.distance = NA
      x$d$geodesic.distance[conns] = real.lengths
    }
    x$d$geodesic.distance.norm = x$d$geodesic.distance/max(x$d$geodesic.distance, na.rm = TRUE)
    relevant.points = subset(x$d, PointNo %in% x$connectors$treenode_id)
    x$connectors$geodesic.distance = relevant.points[match(x$connectors$treenode_id,
                                                           relevant.points$PointNo), ]$geodesic.distance
    x$connectors$geodesic.distance.norm = relevant.points[match(x$connectors$treenode_id,
                                                           relevant.points$PointNo), ]$geodesic.distance.norm
    x
}

#' @export
#' @rdname geodesic_distance
geodesic_distance.neuronlist <- function (x, graph.distance=TRUE, Label = 1, ...) {
  nat::nlapply(x, geodesic_distance.neuron, graph.distance = graph.distance, ...)
}
