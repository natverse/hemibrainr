##########################################################################
################################ meshes ##################################
##########################################################################

#' Get all the region of interest meshes in the hemibrain project
#'
#' @description Get all the ROI meshes as named \code{rgl::mesh3d} objects.
#'
#' @inheritParams hemibrain_read_neurons
#' @inheritParams neuprintr::neuprint_ROIs
#' @inheritParams nat::nlapply
#' @param rois The ROIs to fetch (default \code{NULL} implies all ROIs)
#' @param ... Additional arguments passed to \code{\link{nlapply}}
#' @return a list of \code{rgl::mesh3d} objects
#'
#' @examples
#' \donttest{
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
#' }}
#' @export
#' @seealso \code{\link{hemibrain_skeleton_check}}, \code{\link{hemibrain.surf}}
#' @importFrom neuprintr neuprint_ROIs neuprint_ROI_mesh
hemibrain_roi_meshes <- function(rois=NULL, microns = FALSE, superLevel=NULL, OmitFailures=TRUE, conn=NULL, dataset=NULL, ...){
  if(isTRUE(is.null(rois)))
    rois = unique(neuprint_ROIs(superLevel = superLevel, conn=conn, dataset = dataset))

  fakelist=as.list(structure(rois, .Names=rois))
  hemibrain.rois <- nlapply(fakelist, function(roi) {
    mesh = tryCatch(neuprint_ROI_mesh(roi), error = function(e) NULL)
    if(is.null(mesh)) next
    if(microns){
      mesh*(8/1000)
      nat.templatebrains::regtemplate(mesh) = "JRCFIB2018F"
    }else{
      nat.templatebrains::regtemplate(mesh) = "JRCFIB2018Fraw"
    }
    mesh
    }, OmitFailures=OmitFailures, ...)
  hemibrain.rois
}
