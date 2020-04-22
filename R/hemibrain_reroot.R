######################################################################################
################################ Preprocess neurons ##################################
######################################################################################

#' Re-root a neuron/neurons
#'
#' @description Re-root neurons by predicting their soma location. This works
#' by finding the longest path outside of a given neuropil mesh/meshes, as in insect
#' neurons the soma should be located in a cortex outside of the neuropil proper.
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param meshes a list/a single object of class \code{mesh3d} or \code{hxsurf}.
#' @param ... methods sent to \code{nat::nlapply}
#'
#' @return a \code{nat::neuronlist} or \code{nat::neuron} object
#'
#' @examples
#' \donttest{
#' # Choose some known trouble makers
#' bad.soma = c("5813015982","885788485","5813013913")
#'
#' # Read in these neurons
#' neurons.bs = neuprintr::neuprint_read_neurons(bad.soma)
#'
#' # Re-root
#' neurons = hemibrain_reroot(neurons.bs, meshes = hemibrain.surf)
#'
#' \dontrun{
#' # Let's check that this worked
#' nat::nopen3d()
#' for(n in neurons){
#'      clear3d();
#'      message(n$bodyid)
#'      plot3d(n,col="brown",lwd=2)
#'      points3d(nat::xyzmatrix(n$connectors),col="black")
#'      spheres3d(nat::xyzmatrix(n)[nat::rootpoints(n),],radius=500, alpha=0.5)
#'      p = readline("Done? ")
#'}
#' }}
#' @export
#' @seealso \code{\link{flow_centrality}}
hemibrain_reroot <-function(x, meshes = hemibrain_roi_meshes(), ...)
  UseMethod("hemibrain_reroot")

#' @export
hemibrain_reroot.neuron <- function(x, meshes = hemibrain_roi_meshes(), ...){
  # Find out of volume points
  x$d$roi = NA
  root.orig = nat::rootpoints(x)
  if(!is.hxsurf(meshes)){
    for(roi in names(meshes)){
      inside = nat::pointsinside(nat::xyzmatrix(x$d), surf = meshes[[roi]])
      x$d$roi[inside] = roi
    }
  }else{
    inside = nat::pointsinside(nat::xyzmatrix(x$d), surf = meshes)
    x$d$roi[inside] = deparse(substitute(meshes))
  }
  # Find longest path outside
  outside = rownames(x$d)[is.na(x$d$roi)]
  if(length(outside)==0|is.null(outside)){
    warning("Neuron ", x$bodyid," has no nodes outside of mesh(es), dropping ...")
    x
  }else{
    segs.out = sapply(x$SegList, function(s) sum(s%in%outside))
    longest.out = x$SegList[[which.max(segs.out)]]
    # How many leaf nodes down and upstream of this segment?
    n = nat::as.ngraph(x)
    leaves = nat::endpoints(x)
    downs = igraph::distances(n, v = longest.out, to = leaves, mode = "out")
    downs[is.infinite(downs)] = 0
    downs = apply(downs,1,function(r) sum(r>0))
    down = downs[which.max(downs)]
    ups = igraph::distances(n, v = longest.out, to = leaves, mode = "in")
    ups[is.infinite(ups)] = 0
    ups = apply(ups,1,function(r) sum(r>0))
    up = ups[which.max(ups)]
    if(up<down){
      root = longest.out[length(longest.out)]
    }else{
      root = longest.out[1]
    }
    somid = x$d$PointNo[match(root, 1:nrow(x$d))]
    y = nat::as.neuron(nat::as.ngraph(x$d), origin = somid)
    y = hemibrain_carryover_labels(x=x,y=y)
    y = hemibrain_carryover_tags(x=x,y=y)
    y$connectors = x$connectors
    y$connectors$treenode_id = y$d$PointNo[match(x$connectors$treenode_id, y$d$PointNo)]
    y$soma = y$tags$soma = somid
    y$tags$soma.edit = ifelse(root.orig==root,FALSE,"estimated")
    y = hemibrain_neuron_class(y)
    y
  }
}


#' @export
hemibrain_reroot.neuronlist <- function(x, meshes = hemibrain_roi_meshes(), ...){
  neurons = suppressWarnings(nat::nlapply(x, hemibrain_reroot.neuron, meshes, ...))
  if(sum(!names(x)%in%names(neurons))>0){
    missed = setdiff(names(x),names(neurons))
    warning("Some neurons could not be rerooted, probably do not have points outside meshes: ", paste(missed, collapse = ", "))
  }
  neurons
}


#' Remove incorrectly placed synapses
#'
#' @description With the hemibrain projects somas and their tracts were not
#' separately and explicitly detected. Some false synapse assignment has happened,
#' adding synapses to soma and primary neurite locations, that should not exist
#' based on an examination of the greyscale EM data. This function removes 'out-of-mesh'
#' synapses, you might want to use all ROIs or the hemibrain mesh to remove these 'bad'
#' synapses. In addition, synapses are removed if they are too near the soma / along the
#' primary neurite.
#' @inheritParams flow_centrality
#' @param meshes a list/a single object of class \code{mesh3d} or \code{hxsurf}.
#' @param soma logical, if TRUE it is assumed that the neuron being given has the soma
#' as its route and an intact primary neurite tract. Synapses will then be pruned from both.
#' @param min.nodes.from.soma the minimum number of nodes (geodesic distance) a
#'  synapse can be from the soma. Synapses closes than this will be removed. For
#'  comparable results across neurons, recommended to use \code{nat::resample} on
#'  your neurons before using.
#' @param min.nodes.from.pnt the minimum number of nodes (geodesic distance) a
#'  synapse can be from a point along the primary neurite.
#'  Synapses closes than this will be removed.
#' @param ... methods sent to \code{nat::nlapply}
#'
#' @return a \code{nat::neuronlist} or \code{nat::neuron} object
#'
#' @examples
#' \donttest{
#' # Read in a problematic neuron
#' neuron = neuprint_read_neurons(5813020793)
#'
#' \dontrun{
#' library(nat)
#' plot3d(neuron)
#' points3d(xyzmatrix(neuron[[1]]$connectors))
#' }
#'
#' # Re-root
#' neuron.fixed = hemibrain_remove_bad_synapses(neuron)
#'
#' \dontrun{
#' library(nat)
#' # Let's check that this worked
#' plot3d(neuron.fixed)
#' points3d(xyzmatrix(neuron.fixed[[1]]$connectors), size = 5, col = "green")
#'
#' }}
#' @export
#' @seealso \code{\link{hemibrain_reroot}}
hemibrain_remove_bad_synapses <- function(x,
                                          meshes = hemibrain_roi_meshes(),
                                          soma = TRUE,
                                          min.nodes.from.soma = 100,
                                          min.nodes.from.pnt = 5,
                                          primary.branchpoint = 0.25,
                                          ...) UseMethod("hemibrain_remove_bad_synapses")


#' @export
hemibrain_remove_bad_synapses.neuron <- function(x,
                                                 meshes = hemibrain_roi_meshes(),
                                                 soma = TRUE,
                                                 min.nodes.from.soma = 100,
                                                 min.nodes.from.pnt = 5,
                                                 primary.branchpoint = 0.25,
                                                 ...){
  x.safe = x
  if(!is.null(meshes)){
    x$inside = NA
    if(is.hxsurf(meshes)){
      inside = nat::pointsinside(nat::xyzmatrix(x$connectors), surf = meshes)
      x$connectors$inside[inside] = deparse(substitute(meshes))
    }else{
      for(roi in names(meshes)){
        inside = nat::pointsinside(nat::xyzmatrix(x$connectors), surf = meshes[[roi]])
        x$connectors$inside[inside] = roi
      }
    }
    if(sum(is.na(x$connectors$inside))>0){
      x$connectors = x$connectors[!is.na(x$connectors$inside),]
    }
  }
  if(soma){
    primary.branch.point = primary_branchpoint(x, primary_neurite = TRUE, first = primary.branchpoint)
    pnt = suppressWarnings(unique(unlist(igraph::shortest_paths(nat::as.ngraph(x), nat::rootpoints(x), to = primary.branch.point, mode = "all")$vpath)))
    x$connectors = x$connectors[!x$connectors$treenode_id%in%pnt,]
    syns = unique(x$connectors$treenode_id)
    syns = (1:nrow(x$d))[match(syns,x$d$PointNo)]
    # not within radius of soma
    dists = igraph::distances(nat::as.ngraph(x), v = nat::rootpoints(x), to = syns, mode = "all")
    names(dists) = syns
    x$connectors = x$connectors[x$connectors$treenode_id %in% names(dists)[dists>min.nodes.from.soma],]
    # not within radius of pnt
    dists = igraph::distances(nat::as.ngraph(x), v = pnt, to = syns, mode = "all")
    dists = apply(dists, 2, function(x) sum(x<min.nodes.from.pnt))
    names(dists) = syns
    x$connectors = x$connectors[x$connectors$treenode_id %in% names(dists)[dists==0],]
  }
  removed = nrow(x$connectors) - nrow(x.safe$connectors)
  if(removed>0){
    x$tags$synapses.removed = removed
  }else{
    x$tags$synapses.removed = 0
  }
  x = hemibrain_neuron_class(x)
  x
}

#' @export
hemibrain_remove_bad_synapses.neuronlist <- function(x,
                                                     meshes = hemibrain_roi_meshes(),
                                                     soma = TRUE,
                                                     min.nodes.from.soma = 125,
                                                     min.nodes.from.pnt = 5,
                                                     primary.branchpoint = 0.25,
                                                     ...){
  nat::nlapply(x,
                        hemibrain_remove_bad_synapses.neuron,
                        meshes = meshes,
                        soma = soma,
                        min.nodes.from.soma=min.nodes.from.soma,
                        min.nodes.from.pnt=min.nodes.from.pnt,
                        primary.branchpoint = primary.branchpoint,
                        ...)
}


#' Pre-process hemibrain skeletons
#'
#' @description Call \code{\link{hemibrain_reroot}} on neurons read from
#' the hemibrain neuprintr project, that have no marked soma, and
#' \code{\link{hemibrain_remove_bad_synapses}} on all neurons.
#' @inheritParams hemibrain_remove_bad_synapses
#' @inheritParams hemibrain_reroot
#' @param OmitFailures Whether to omit neurons for which FUN gives an error.
#' The default value (NA) will result in nlapply stopping with an error
#' message the moment there is an error.
#'
#' @return a \code{nat::neuronlist} or \code{nat::neuron} object
#'
#' @examples
#' \donttest{
#' # Choose some known trouble makers
#' ids = c("5813056323", "579912201", "5813015982", "973765182", "885788485",
#' "915451074", "5813032740", "1006854683", "5813013913", "5813020138",
#' "853726809", "916828438", "5813078494", "420956527", "486116439",
#' "573329873", "5813010494", "5813040095", "514396940", "665747387",
#' "793702856", "451644891", "482002701", "391631218", "390948259",
#' "390948580", "452677169", "511262901", "422311625", "451987038"
#' )
#'
#' # just use a subset for documentation purposes
#' ids=ids[1:5]
#'
#' # Read in these neurons
#' neurons = neuprintr::neuprint_read_neurons(ids)
#'
#' # Get all the roi meshes
#' hemibrain.rois = hemibrain_roi_meshes()
#' ## Using this as the argument for 'meshes' will also
#' ## give you the ROI for each point and synapse in a neuron's skeleton!!
#'
#' # Re-root
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.rois)
#'
#' \dontrun{
#' # Let's check that this worked
#' nat::nopen3d()
#' for(n in neurons.checked){
#'      clear3d();
#'      message(n$bodyid)
#'      plot3d(n,col="brown",lwd=2)
#'      points3d(nat::xyzmatrix(n$connectors),col="black")
#'      spheres3d(nat::xyzmatrix(n)[nat::rootpoints(n),],radius=500, alpha=0.5)
#'      p = readline("Done? ")
#'}
#' }}
#' @export
#' @seealso \code{\link{hemibrain_reroot}}
#' @importFrom nat as.neuronlist
hemibrain_skeleton_check <- function(x, # as read by neuprint_read_neurons
                                     meshes = hemibrain_roi_meshes(),
                                     min.nodes.from.soma = 100,
                                     min.nodes.from.pnt = 5,
                                     OmitFailures = FALSE,
                                     ...){

  # Re-root somas where necessary
  x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  x.nosoma = x[!x[,"soma"]]
  x.soma = nat::setdiff(x,x.nosoma)
  x.soma = hemibrain_settags(x.soma, soma.edit = rep(FALSE,length(x.soma)))
  if(length(x.nosoma)){
    message("Re-rooting ", length(x.nosoma), " neurons without a soma")
    x.estsoma = hemibrain_reroot(x = x.nosoma, meshes = meshes, OmitFailures = OmitFailures, ...)
    nams = intersect(names(x),c(names(x.soma),names(x.estsoma)))
    x.new = c(x.soma, x.estsoma)[nams]
  } else {
    x.new = x.soma
  }

  # Remove erroneous synapses, out of mesh and on pnt/soma
  message("Removing synapses at somas and along primary neurite for ", length(x.new), " neurons")
  x.goodsyn = hemibrain_remove_bad_synapses(x.new,
                                            meshes = meshes,
                                            soma = TRUE,
                                            min.nodes.from.soma = min.nodes.from.soma,
                                            min.nodes.from.pnt = min.nodes.from.pnt,
                                            OmitFailures = OmitFailures,
                                            ...)

  # Add new info to meta-data
  x.goodsyn = metadata_add_tags(x.goodsyn)
  x.goodsyn
}
