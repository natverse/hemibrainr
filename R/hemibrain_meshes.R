##########################################################################
################################ meshes ##################################
##########################################################################

#' Get some or all the region of interest meshes in the hemibrain project
#'
#' @description This function can retrieve some or all the hemibrain project ROI
#'   meshes as \code{rgl::\link{mesh3d}} objects. They are returned in a named
#'   \code{\link{neuronlist}} object. This allows further subsets and colours to
#'   be applied via the standard \code{nat::\link{plot3d.neuronlist}} semantics.
#'   See the examples for further details.
#'
#' @inheritParams hemibrain_read_neurons
#' @inheritParams neuprintr::neuprint_ROIs
#' @inheritParams nat::nlapply
#' @param rois The ROIs to fetch (default \code{NULL} implies all ROIs)
#' @param ... Additional arguments passed to \code{\link{nlapply}}
#' @return a \code{neuronlist} of \code{rgl::mesh3d} objects
#'
#' @examples
#' \dontrun{
#' # Get all the roi meshes
#' hemibrain.rois = hemibrain_roi_meshes()
#'
#' # Read in these neurons
#' cols = rainbow(length(hemibrain.rois))
#' for(m in 1:length(hemibrain.rois)){
#'     plot3d(hemibrain.rois[[m]],alpha = 0.5, col = cols[m])
#' }
#'
#' # Alternatively plot a subset all at once with rainbow colours
#' # making use of nat::plot3d.neuronlist
#' plot3d(hemibrain.rois, 1:3, add=TRUE)
#'
#' # If you have neuprintr >= 1.0 you can easily find groups of meshes
#' # find just the terminal nodes of the MB ROI hierarchy
#' library(nat)
#' mbterms=nat::endpoints(neuprint_ROI_hierarchy(root='MB(R)', rval='graph'))
#' mb.rois=hemibrain_roi_meshes(mbterms)
#' plot3d(mb.rois, add=TRUE)
#' }
#' @export
#' @seealso \code{\link{hemibrain_skeleton_check}}, \code{\link{hemibrain.surf}},
#' \code{neuprintr::\link[neuprintr]{neuprint_ROI_mesh}},
#' \code{neuprintr::\link[neuprintr]{neuprint_ROI_hierarchy}}
#' @importFrom neuprintr neuprint_ROIs neuprint_ROI_mesh
hemibrain_roi_meshes <- function(rois=NULL, microns = FALSE, superLevel=NULL, OmitFailures=TRUE, conn=NULL, dataset=NULL, ...){
  if(isTRUE(is.null(rois)))
    rois = unique(neuprint_ROIs(superLevel = superLevel, conn=conn, dataset = dataset))

  fakelist=as.list(structure(rois, .Names=rois))
  hemibrain.rois <- nlapply(fakelist, function(roi) {
    mesh = neuprint_ROI_mesh(roi)
    if(microns){
      mesh <- mesh*(8/1000)
      nat.templatebrains::regtemplate(mesh) = "JRCFIB2018F"
    }else{
      nat.templatebrains::regtemplate(mesh) = "JRCFIB2018Fraw"
    }
    mesh
    }, OmitFailures=OmitFailures, ...)
  hemibrain.rois
}
