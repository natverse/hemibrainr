#' Read neurons from the hemibrain connectome project
#'
#' @description Read neurons from the hemibrain connectome project. This function
#' uses the package \code{neuprintr} to read neurons. It then uses
#' \code{\link{hemibrain_flow_centrality}} to re-root and split neurons into putative
#' axons and dendrites. Optionally, it may also convert neurons from their raw voxel
#' space (as they are stored in NeuPrint) to microns.
#'
#' @param x a vector of bodyids that can be read from \url{'https://neuprint.janelia.org/'}.
#' @param microns convert dimensions from raw voxels into microns (template brain: \code{JRCFIB2018F}, else \code{JRCFIB2018Fraw}).
#' @param ... methods passed to \code{neuprintr::neuprint_read_neurons}, \code{\link{hemibrain_remove_bad_synapses}}
#'  and \code{\link{hemibrain_flow_centrality}}
#'
#'
#' @inherit flow_centrality return
#'
#' @examples
#' \donttest{
#'
#' # Choose neurons
#' ## In this case some antennal lobe local neurons
#' al.local.neurons = c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#' "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#' "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#' "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#' "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#' "1794037618", "5813018729", "2013333009")
#'
#' # Get neurons
#' neurons = hemibrain_read_neurons(al.local.neurons)
#'
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons, WithConnectors = TRUE)
#' }}
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}},\code{\link{hemibrain_remove_bad_synapses}}
hemibrain_read_neurons<-function(x, microns = TRUE, ...){
  neurons = neuprintr::neuprint_read_neurons(x, ...)
  neurons.flow = hemibrain_flow_centrality(neurons, ...)
  neurons.flow = hemibrain_remove_bad_synapses(neurons.flow, ...)
  hemibrain_metrics = hemibrain_metrics[,!colnames(hemibrain_metrics)%in%colnames(neurons.flow[,])]
  df = cbind(neurons.flow[,], hemibrain_metrics[names(neurons.flow),])
  rownames(df) = names(neurons.flow)
  if(microns){
    neurons.flow = scale_neurons(neurons.flow, scaling = (8/1000))
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018F"
  }else{
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018Fraw"
  }
  neurons.flow[,] = df
  neurons.flow
}

# hidden
scale_neurons <-function(x, scaling = (8/1000), ...) UseMethod("scale_neurons")
scale_neurons.neuron <- function(x, scaling, ...){
  nat::xyzmatrix(x$d) = nat::xyzmatrix(x$d)*scaling
  nat::xyzmatrix(x$connectors) = nat::xyzmatrix(x$connectors)*scaling
  x
}
scale_neurons.neuronlist = function(x, scaling = (8/1000), ...){
  nat::nlapply(x,scale.neuron, scaling = scaling, ...)
}

