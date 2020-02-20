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




#' Axon-dendrite split points calculated for the whole dataset
#'
#' Axon-dendrite split points calculated for the whole hemibrain dataset,
#' hemibrain:v1.0.1. The fiel used to generate this can be found in this package,
#' as \code{data-raw/hemibrain_split_all_neurons.R}. Because synapses can be unevently
#' distributed in this dataset across brain regions, after compartment detection,
#' the dendrite was chosen as the neurite closest to a neuron's primary branch point.
#'
#' @references Bates, A. S., P. Schlegel, R. J. V. Roberts, N. Drummond, I. F. M. Tamimi, R. Turnbull, X. Zhao, et al. 2020.
#' “Complete Connectomic Reconstruction of Olfactory Projection Neurons in the Fly Brain.”
#' bioRxiv. https://doi.org/10.1101/2020.01.19.911453.
#'
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.19.911453v1}
#' @seealso \code{\link{hemibrain_flow_centrality}}, \code{\link{hemibrain_splitpoints}}, \code{\link{flow_centrality}}, \code{\link{hemibrain_use_splitpoints}}
#' @examples
#' \donttest{
#' View(hemibrain_splitpoints_distance)
#' }
"hemibrain_splitpoints_distance"
