##########################################################################
################################ meshes ##################################
##########################################################################

#' Get all the region of interest meshes in the hemibrain project
#'
#' @description Get all the ROI meshes as named \code{rgl::mesh3d} objects.
#'
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
#' @seealso \code{\link{hemibrain_skeleton_check}}
hemibrain_roi_meshes <- function(){
  rois = neuprintr::neuprint_ROIs(superLevel = NULL)
  hemibrain.rois = list()
  for(roi in rois){
    mesh = tryCatch(neuprintr::neuprint_ROI_mesh(roi), error = function(e) NULL)
    hemibrain.rois[[roi]] = mesh
  }
  hemibrain.rois
}

