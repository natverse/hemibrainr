######################################################################################
################################ Preprocess neurons ##################################
######################################################################################

#' Re-root (a) hemibrain/flywire neuron/neurons
#'
#' @description Re-root neurons by predicting their soma location. Somas have been manually tagged by the FlyEM project.
#' However, in some cases the roots are wrong, or somas are outside of the volume. The Fly Connectome team at the University of Cambridge has
#' manually tagged better root points for these neurons. Soma locations, updated by this work, as saved as \code{\link{hemibrain_somas}} in this package.
#' A user can read somas from here, or get a bleeding edge version from
#' the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}. Alternatively, they
#' can 'estimated' a good root-point. This works by finding the longest path outside of a given neuropil mesh/meshes, as in insect
#' neurons the soma should be located in a cortex outside of the neuropil proper.
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param method whether to use the manually curated soma list or estimate soma location.
#' @param googlesheet logical, whether to read soma locations from the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}
#' if \code{method == "manual"}.
#' @param hemibrain_somas a \code{data.frame} that gives soma locations for hemibrain neurons. See the default, \code{\link{hemibrain_somas}}. If \code{googelsheet} is \code{TRUE} this is read fresh
#' from the hemibrain Google team drive overseen by the Drosophila Connectomics group.
#' @param flywire_nuclei a \code{data.frame} that gives auto-detected nuclei locations in the Flywire dataset. This is provided by
#' @param meshes a list/a single object of class \code{mesh3d} or \code{hxsurf}. Only used for estimation.
#' If \code{NULL} then \code{hemibrain_roi_meshes} is called.
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
#' # Let's check that this worked. Old root in red, new in green
#' nat::nopen3d()
#' for(i in 1:length(neurons)){
#'      n = neurons[[i]]
#'      bs = neurons.bs[[i]]
#'      clear3d();
#'      message(n$bodyid)
#'      plot3d(n,col="brown",lwd=2)
#'      points3d(nat::xyzmatrix(n$connectors),col="black")
#'      spheres3d(nat::xyzmatrix(n)[nat::rootpoints(n),],radius=500, alpha=0.5, col = "green")
#'      spheres3d(nat::xyzmatrix(bs)[nat::rootpoints(bs),],radius=500, alpha=0.5, col = "red")
#'      p = readline("Done? ")
#'}
#' }}
#' @export
#' @seealso \code{\link{flow_centrality}}, \code{\link{hemibrain_somas}}
hemibrain_reroot <-function(x,
                            method = c("manual","estimated"),
                            meshes = NULL,
                            googlesheet = FALSE,
                            hemibrain_somas = hemibrainr::hemibrain_somas,
                            ...) UseMethod("hemibrain_reroot")

#' @export
hemibrain_reroot.neuron <- function(x,
                                    method = c("manual","estimated"),
                                    meshes = NULL,
                                    googlesheet = FALSE,
                                    hemibrain_somas = hemibrainr::hemibrain_somas,
                                    ...){
  method = match.arg(method)
  # Find out of volume points
  root.orig = nat::rootpoints(x)
  if(method == "manual"){
    if(googlesheet){
      selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E"
      hemibrain_somas = googlesheets4::read_sheet(ss = selected_file, sheet = "somas")
      hemibrain_somas = as.data.frame(hemibrain_somas, stringsAsFactors = FALSE)
      rownames(hemibrain_somas) = hemibrain_somas$bodyid
    }
    bis = as.character(x$bodyid)
    if(is.null(bis)){
      stop("Neuron must have its bodyid at x$bodyid. If this is not the case, try running x = add_field_seq(x,x[,'bodyid'],field='bodyid')")
    }else if (!bis%in%rownames(hemibrain_somas)){
      warning("Neuron ",bis," not in saved soma locations. Estimating soma")
      root = root.orig
    }else{
      som = hemibrain_somas[bis,]
      root = nabor::knn(query = nat::xyzmatrix(som), data = nat::xyzmatrix(x$d), k = 1)$nn.idx
    }
    somid = x$d$PointNo[match(root, 1:nrow(x$d))]
    y = nat::as.neuron(nat::as.ngraph(x$d), origin = somid)
    y$tags$soma.edit = ifelse(root.orig==root,FALSE,"manual")
  }
  if (method == "estimated"){
    if(is.null(meshes)){
      meshes = hemibrain_roi_meshes()
    }
    x$d$roi = NA
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
      y$tags$soma.edit = ifelse(root.orig==root,FALSE,"estimated")
    }
  }

  # Return
  y = hemibrain_carryover_labels(x=x,y=y)
  y = hemibrain_carryover_tags(x=x,y=y)
  y$connectors = x$connectors
  y$connectors$treenode_id = y$d$PointNo[match(x$connectors$treenode_id, y$d$PointNo)]
  y$soma = y$tags$soma = somid
  y = hemibrain_neuron_class(y)
  y
}

#' @export
hemibrain_reroot.neuronlist <- function(x,
                                        method = c("manual","estimated"),
                                        meshes = NULL,
                                        googlesheet = FALSE,
                                        ...){
  method = match.arg(method)
  if(googlesheet){
    selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E"
    hemibrain_somas = googlesheets4::read_sheet(ss = selected_file, sheet = "somas")
    hemibrain_somas = as.data.frame(hemibrain_somas, stringsAsFactors = FALSE)
    rownames(hemibrain_somas) = hemibrain_somas$bodyid
  }
  if(method == "estimated"){
    if(is.null(meshes)){
      meshes = hemibrain_roi_meshes()
    }
  }
  x = tryCatch(add_field_seq(x,x[,"bodyid"],field="bodyid"),
               error = function(e) add_field_seq(x,names(x),field="bodyid"))
  neurons = suppressWarnings(nat::nlapply(x, hemibrain_reroot.neuron, meshes = meshes, hemibrain_somas = hemibrain_somas, googlesheet = FALSE, ...))
  if(sum(!names(x)%in%names(neurons))>0){
    missed = setdiff(names(x),names(neurons))
    warning("Some neurons could not be rerooted, probably do not have points outside meshes: ", paste(missed, collapse = ", "))
  }
  neurons
}

#' @export
flywire_reroot <- function(x,
                           flywire_nuclei = fafbseg::flywire_nuclei(),
                           ...) UseMethod("flywire_reroot")

#' @export
flywire_reroot.neuron <- function(x, flywire_nuclei = fafbseg::flywire_nuclei(), ...){
  root.id = as.character(x$root_id)
  flywire_nuclei$pt_root_id = as.character(flywire_nuclei$pt_root_id)
  if(is.null(root.id)){
    stop("no root_id at x$rootid for given neuron")
  }
  if(!root.id%in%flywire_nuclei$pt_root_id){
    warning(root.id, " not in flywire_nuclei")
    y = x
  }else{
    flywire.nucleus = subset(flywire_nuclei, flywire_nuclei$pt_root_id == root.id)[1,]
    som = matrix(flywire.nucleus$pt_position[[1]], ncol = 3)
    root = nabor::knn(query = som, data = nat::xyzmatrix(x$d), k = 1)$nn.idx
    somid = x$d$PointNo[match(root, 1:nrow(x$d))]
    y = nat::as.neuron(nat::as.ngraph(x$d), origin = somid)
    y = hemibrain_carryover_labels(x=x,y=y)
    y = hemibrain_carryover_tags(x=x,y=y)
    y$connectors = x$connectors
    y$connectors$treenode_id = y$d$PointNo[match(x$connectors$treenode_id, y$d$PointNo)]
    flywire.nucleus$treenode_id = somid
    y$soma = flywire.nucleus
    if(!is.null(y$d$Label)){
      y$d$Label[somid] = 1
    }
  }
  # return
  y = hemibrain_neuron_class(y)
  y
}

#' @export
flywire_reroot.neuronlist <- function(x, flywire_nuclei = fafbseg::flywire_nuclei(), ...){
  x = tryCatch(add_field_seq(x,x[,"root_id"],field="root_id"),
               error = function(e) add_field_seq(x,names(x),field="root_id"))
  y = nat::nlapply(X = x, FUN = flywire_reroot.neuron, flywire_nuclei = flywire_nuclei, ...)
  y
}

#' Remove incorrectly placed synapses
#'
#' @description With the hemibrain projects somas and their tracts were not
#'   separately and explicitly detected. Some false synapse assignment has
#'   happened, adding synapses to soma and primary neurite locations, that
#'   should not exist based on an examination of the greyscale EM data. This
#'   function removes 'out-of-mesh' synapses, you might want to use all ROIs or
#'   the hemibrain mesh to remove these 'bad' synapses. In addition, synapses
#'   are removed if they are too near the soma / along the primary neurite.
#' @inheritParams flow_centrality
#' @param meshes a list/a single object of class \code{mesh3d} or \code{hxsurf}.
#'   Defaults to \code{\link{hemibrain.surf}}. See examples for alternatives.
#' @param soma logical, if TRUE it is assumed that the neuron being given has
#'   the soma as its route and an intact primary neurite tract. Synapses will
#'   then be pruned from both.
#' @param min.nodes.from.soma the minimum number of nodes (geodesic distance) a
#'   synapse can be from the soma. Synapses closes than this will be removed.
#'   For comparable results across neurons, recommended to use
#'   \code{nat::resample} on your neurons before using.
#' @param min.nodes.from.pnt the minimum number of nodes (geodesic distance) a
#'   synapse can be from a point along the primary neurite. Synapses closes than
#'   this will be removed.
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
#' # remove any synapses outside of this surface (or on the cell body fibre)
#' neuron.fixed = remove_bad_synapses(neuron, meshes=hemibrain.surf)
#' \dontrun{
#' # previously we used all the individual hemibrain meshes. This may produce
#' # slightly different results but using *hemibrain.surf* is much faster
#' neuron.fixed = remove_bad_synapses(neuron, meshes=hemibrain_roi_meshes())
#' }
#' \dontrun{
#' library(nat)
#' # Let's check that this worked
#' plot3d(neuron.fixed)
#' points3d(xyzmatrix(neuron.fixed[[1]]$connectors), size = 5, col = "green")
#'
#' }
#' \dontrun{
#' # Now an example with flywire neurons.

#' }
#'}
#' @export
#' @seealso \code{\link{hemibrain_reroot}}
remove_bad_synapses <- function(x,
                                meshes = NULL, # hemibrainr::hemibrain.surf,
                                soma = TRUE,
                                min.nodes.from.soma = 100,
                                min.nodes.from.pnt = 5,
                                primary.branchpoint = 0.25,
                                ...) UseMethod("remove_bad_synapses")



# # Get our exemplar FAFB neurons
# nx <- xform_brain(elmr::dense_core_neurons, ref="FlyWire", sample="FAFB14")
#
# # Get flywire IDs
# choose_segmentation("flywire")
# xyz <- xyzmatrix(nx)
# ids <- unique(flywire_xyz2id(xyz[sample(1:nrow(xyz),100),]))
#
# # Skeletonise
# neurons <- fafbseg::skeletor(ids,
#                    method = "wavefront",
#                    save.obj = NULL,
#                    mesh3d = FALSE,
#                    waves = 1)
#
# # Re-root the soma
# # Get flywire nuclei
# flywire_n <- fafbseg::flywire_nuclei()
# flywire_n <- as.data.frame(flywire_n)
# flywire_n$pt_root_id <- as.character(flywire_n$pt_root_id)
# neurons.rerooted <- flywire_reroot(neurons, .parallel = FALSE, flywire_nuclei = flywire_n)
#
# # Add synapses
# neurons.syn <- fafbseg::flywire_neurons_add_synapses(x = neurons.rerooted,
#                                                      cleft.threshold = 50,
#                                                      transmitters = TRUE,
#                                                      local = NULL,
#                                                      Verbose = FALSE,
#                                                      OmitFailures = FALSE,
#                                                      .parallel = FALSE)
#
# # Split into axon and dendrite
# neurons.syn.split <- flow_centrality(neurons.syn)
#
# # Flag bad synapses
# neurons.syn.split.flagged <- remove_bad_synapses(neurons.syn.split,
#                                                  soma = TRUE,
#                                                  min.nodes.from.soma = 100,
#                                                  min.nodes.from.pnt = 5,
#                                                  primary.branchpoint = 0.25)
#
# # NT predictions in good versus bad synapses
# syns <- nlapply(neurons.syn.split.flagged, function(x) x$connectors)
# syns.df <- do.call(rbind, syns)
# table(syns.df$Label, syns.df$status)

#' @export
remove_bad_synapses.neuron <- function(x,
                                       meshes = NULL, # hemibrainr::hemibrain.surf,
                                       soma = TRUE,
                                       min.nodes.from.soma = 50,
                                       min.nodes.from.pnt = 5,
                                       primary.branchpoint = 0.25,
                                       ...){
  x.safe = x
  split = !is.null(x$AD.segregation.index)
  if(is.null(x$connectors$status)){
    x$connectors$status='good'
  }
  if(!is.null(meshes)){
    x$connectors$inside = NA
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
  if(split){
    soma.nodes <- x$connectors$Label %in% c("soma",'1',1)
    pnt.nodes <- x$connectors$Label %in% c("pnt",'primary_neurite','primary.neurite','7',7)
    pd.nodes <- x$connectors$Label %in% c("pnt",'primary_dendrite','primary.dendrite','4',4)
    x$connectors[soma.nodes,"status"]="on_soma"
    x$connectors[pnt.nodes,"status"]="on_pnt"
    x$connectors[pd.nodes,"status"]="on_pd"
  }else{
    primary.branch.point = primary_branchpoint(x, primary_neurite = TRUE, first = primary.branchpoint)
    pnt = suppressWarnings(unique(unlist(igraph::shortest_paths(graph = nat::as.ngraph(x),
                                                                from = nat::rootpoints(x),
                                                                to = primary.branch.point,
                                                                mode = "all")$vpath)))
    x$connectors[x$connectors$treenode_id%in%pnt,"status"]="near_pnt"
    syns = unique(x$connectors$treenode_id)
    syns = (1:nrow(x$d))[match(syns,x$d$PointNo)]
    # not within radius of pnt
    dists = igraph::distances(nat::as.ngraph(x), v = pnt, to = syns, mode = "all")
    dists = apply(dists, 2, function(x) sum(x<min.nodes.from.pnt))
    names(dists) = syns
    x$connectors[x$connectors$treenode_id %in% names(dists)[as.numeric(dists)==0],"status"]="on_pnt"
  }
  if(soma){
    # not within radius of soma
    syns = unique(x$connectors$treenode_id)
    syns = (1:nrow(x$d))[match(syns,x$d$PointNo)]
    dists = igraph::distances(graph = nat::as.ngraph(x),
                              v = nat::rootpoints(x),
                              to = syns,
                              mode = "all")
    names(dists) = syns
    x$connectors[x$connectors$treenode_id %in% names(dists)[as.numeric(dists)<=min.nodes.from.soma],'status']="near_soma"
  }
  removed = sum(x$connectors$status!="good", na.rm = TRUE)
  if(removed>0){
    x$tags$synapses.removed = removed
    message('flagged ', removed," of ",  nrow(x$connectors), " bad synapses on skeleton, see x$connectors$status")
  }else{
    x$tags$synapses.removed = 0
  }
  x <- hemibrain_neuron_class(x)
  x
}

#' @export
remove_bad_synapses.neuronlist <- function(x,
                                           meshes = NULL, #hemibrainr::hemibrain.surf,
                                           soma = TRUE,
                                           min.nodes.from.soma = 125,
                                           min.nodes.from.pnt = 5,
                                           primary.branchpoint = 0.25,
                                           ...){
  nat::nlapply(x,
              remove_bad_synapses.neuron,
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
#' \code{\link{remove_bad_synapses}} on all neurons. Does not prune synapses
#' from soma positions / along the primary neurite, but does remove synapses outside
#' of the hemibrain volume.
#' @inheritParams remove_bad_synapses
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
#' ids=ids[1:3]
#'
#' # Read in these neurons
#' neurons = neuprintr::neuprint_read_neurons(ids)
#' \dontrun{
#' # Get all the roi meshes
#' hemibrain.rois = hemibrain_roi_meshes()
#' ## Using this as the argument for 'meshes' will also
#' ## give you the ROI for each point and synapse in a neuron's skeleton!!
#'
#' # Re-root
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.rois)
#' }
#'
#' # Alternatively, we can do this faster using a single surface model for the
#' # whole hemibrain
#'
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.surf)
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
                                     meshes = hemibrainr::hemibrain.surf,
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
    x.estsoma = hemibrain_reroot(x = x.nosoma, meshes = meshes, OmitFailures = OmitFailures)
    nams = intersect(names(x),c(names(x.soma),names(x.estsoma)))
    x.new = c(x.soma, x.estsoma)[nams]
  } else {
    x.new = x.soma
  }

  # Remove erroneous synapses, out of mesh and on pnt/soma
  message("Removing synapses outside the hemibrain neuropil volume for ", length(x.new), " neurons")
  x.goodsyn = remove_bad_synapses(x.new,
                                  meshes = meshes,
                                  soma = FALSE,
                                  min.nodes.from.soma = min.nodes.from.soma,
                                  min.nodes.from.pnt = min.nodes.from.pnt,
                                  OmitFailures = OmitFailures)

  # Add new info to meta-data
  x.goodsyn = metadata_add_tags(x.goodsyn)
  x.goodsyn
}

#' Cut neurons to the hemibrain volume
#'
#' @description Cut neurons from another data source to the hemibrain volume.
#'
#' @param x a \code{nat::neuronlist} object.
#' @param brain the brainspace to which \code{x} is registered.
#' @param scale scale \code{x} by a factor of \code{scale}.
#' @param mirror whether or not to mirror the hemibrain volume, to which \code{x} is cut.
#' @param ... methods supplied to \code{nat::nlapply}
#'
#' @return a \code{nat::neuronlist}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' Cell07PNs.cut = hemibrain_cut(Cell07PNs.cut, brain = "FCWB")
#' }}
#' @export
#' @seealso \code{\link{hemibrain_reroot}}
hemibrain_cut <- function(x,
                          scale = 1,
                          brain = "FAFB14",
                          mirror = FALSE,
                          ...){
  fafb.cut.1 = structure(c(312048, 601733, 71265, 319018, 4315, 270859), .Dim = c(2L, 3L))
  fafb.cut.2 = structure(c(519300, 601733, 259570, 319018, 0, 300000), .Dim = c(2L, 3L))
  if(mirror){
    fafb.cut.1 = nat.jrcbrains::mirror_fafb(fafb.cut.1)
    fafb.cut.2 = nat.jrcbrains::mirror_fafb(fafb.cut.2)
  }
  if(brain!="FAFB14"){
    fafb.cut.1 = nat.templatebrains::xform_brain(fafb.cut.1, sample = "FAFB14", reference = brain)
    fafb.cut.2 = nat.templatebrains::xform_brain(fafb.cut.2, sample = "FAFB14", reference = brain)
  }
  y = nat::nlapply(x,
                   subbbx,
                   bbx = fafb.cut.1,
                   scale = scale,
                   ret='inside',
                   ...)
  # Remove a chunk of left antennal lobe that's missing in hemibrain
  z = nat::nlapply(y,
                   subbbx,
                   bbx = fafb.cut.2,
                   scale = scale,
                   ret='outside',
                   ...)
  # drop any FAFB neurons that don't have at least 5 vertices left
  z = z[nat::nvertices(z)>=5]
  z
}

# Cut neurons by bounding box
subbbx <- function(n.dps, bbx, scale, ret='inside'){
  points = nat::xyzmatrix(n.dps)
  bbx = nat::boundingbox(scale(bbx, scale = rep(scale,3), center = FALSE))
  inside = nat::pointsinside(points,bbx)
  # or subset(n.dps, inside, invert=ret != 'inside')
  if (ret == 'inside') {
    subset(n.dps, inside)
  }
  else {
    subset(n.dps,!inside)
  }
}

