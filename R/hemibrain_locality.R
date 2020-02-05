##################################################################################
################################ Locality Score ##################################
##################################################################################

hemibrain_locality.neuron <- function(x, resample = 1000){

  # Extract compartment cable
  axon = axonic_cable(x)
  dend = dendritic_cable(x)

  # Resample compartments
  if(!is.null(resample)){
    axon.res = rnat::resample(axon, stepsize = resample)
    dend.res = rnat::resample(dend, stepsize = resample)
  }else{
    axon.res = axon
    dend.res = dend
  }

  # Calculate compartment overlap scores
  oad = overlap_score(axon.res, dend.res)
  oda = overlap_score(dend.res, axon.res)

  #


}







#' Generate a connectivity matrix based on euclidean distance between points
#'
#' @description Generates an 'overlap matrix' of overlap scores between neurons in the \code{output.neurons} and \code{input.neurons} pools.
#' For every point in a given neuron in \code{output.neurons}, a distance score is calculated to every point in a neuron in \code{input.neurons}.
#' The sum of this score is added to the final output matrix. The score is calculated as \code{e(-d^2/2δ^2)}, where d is the euclidean distance between the two points,
#' and δ is the expected distance in um that is considered 'close'. It is recommended that the user resamples neurons before use, using \code{\link{resample}}.
#'
#' @param output.neurons first set of neurons
#' @param input.neurons second set of neurons
#' @param delta the distance (in um) at which a synapse might occur
#' @param progress whether or not to have a progress bar
#'
#' @examples
#' \dontrun{
#' # Calculate how much some neurons overlap with one another
#' ## Example requires the package nat.flybrains
#' Cell07PNs_overlap = overlap(output.neurons = nat.flybrains::Cell07PNs,
#' input.neurons = nat.flybrains::Cell07PNs)
#'
#' ## Plot the results
#' heatmap(Cell07PNs_overlap)
#' }
#' @return a matrix of overlap scores
#' @seealso \code{\link{hemibrain_locality}}
#' @export
overlap_score <- function(output.neurons, input.neurons, delta =1, progress = TRUE){
  score.matrix = matrix(0,nrow = length(output.neurons),ncol = length(input.neurons))
  rownames(score.matrix) = names(output.neurons)
  colnames(score.matrix) = names(input.neurons)
  for (n in 1:length(output.neurons)){
    a = xyzmatrix(output.neurons[[n]])
    input.neurons.d = nlapply(input.neurons, xyzmatrix, .progress = "none")
    s = sapply(input.neurons.d, function(x)sum(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)))) # Score similar to that in Schlegel et al. 2015
    score.matrix[n,] = s
    if(progress){
      nat_progress(x = n/length(output.neurons)*100, message = "calculating overlap")
    }
  }
  score.matrix
}
