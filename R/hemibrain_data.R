#' The hemibrain mesh
#'
#' A 3D mesh of the hemibrain volume as a
#' \code{\link[nat]{hxsurf}} object. Created by
#' Philipp Schlegel from the individual ROIs in
#' FlyEM's hemibrain project.
#'
#' @source \url{https://neuprint.janelia.org/}
#' @seealso \code{\link{hemibrain_roi_meshes}}
#' @examples
#' \donttest{
#' rgl::wire3d(hemibrain)
#' rgl::shade3d(hemibrain, alpha=.3)
#' rgl::shade3d(hemibrain, alpha=.3, col='grey')
#' }
"hemibrain.surf"
