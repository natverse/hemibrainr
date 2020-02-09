#' The hemibrain mesh
#'
#' A 3D mesh of the hemibrain volume as a \code{\link[nat]{hxsurf}} object.
#' Created by Philipp Schlegel from the union of individual ROIs in FlyEM's
#' hemibrain project.
#'
#' @details note that hemibrain coordinate system does not have the
#'   anterior-posterior axis aligned with the Z axis.
#'
#' @source \url{https://neuprint.janelia.org/}
#' @seealso \code{\link{hemibrain_roi_meshes}}
#' @examples
#' \donttest{
#' library(nat)
#' plot3d(hemibrain.surf, alpha=.3, col='grey')
#'
#' # can also convert to rgl mesh3d object for more control
#' hemibrain <- as.mesh3d(hemibrain.surf)
#' rgl::wire3d(hemibrain)
#' rgl::shade3d(hemibrain, alpha=.3)
#' rgl::shade3d(hemibrain, alpha=.3, col='grey')
#' }
"hemibrain.surf"
