##############################################################################################
################################ Access Data on GoogleDrive ##################################
##############################################################################################

#' Read precomputed information from the hemibrain Googledrive
#'
#' @description Read precomputed data available on the Hemibrain Google Team Drive. This includes all synapses, neuron-neuron connections
#' and an edgelist for all hemibrain neurons, broken down by axon and dendrite assigments.
#'
#' @param savedir directory in which \code{nat::neuronlistfh} object as been saved using \code{\link{hemibrain_download_neurons}}. This acts as a database,
#' from which hemibrain neurons can be read. If \code{TRUE} your default save directory is used, which is stored as: \code{options()$hemibrain_data}
#' @param local logical, whether to try to read locally saved neurons (by default at: \code{options()$hemibrain_data}) or neurons from Google Drive (\code{options()$Gdrive_hemibrain_data}).
#' @param neuron.split read saved neurons spit in which way? Folder names indicative of arguments passed to \code{\link{flow_centrality}}.
#' and \code{\link{hemibrain_flow_centrality}}
#'
#' @return a \code{data.frame} or character vector
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' #  All neurpritn IDs for neurons that have a split precomputed
#' ids = hemibrain_neuron_bodyids()
#'
#' # Connectivity edgelist, broken down by axon/dendrite
#' elist = hemibrain_edgelist()
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}}
#' @rdname hemibrain_googledrive
#' @export
hemibrain_neuron_bodyids <- function(savedir = TRUE,
                                        local = FALSE){
  savedir = good_savedir(savedir = savedir,local = local)
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neuron_bodyids", folder = "hemibrain_neurons/")
  gcsv = read.csv(gfile)
  as.character(gcsv$x)
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_elist <- function(savedir = TRUE, local = FALSE,
                            neuron.split = c("polypre_centrifugal_synapses",
                                             "polypre_centrifugal_distance")){
  neuron.split = match.arg(neuron.split)
  savedir = good_savedir(savedir = savedir,local = local)
  folder = paste0("hemibrain_neurons/",neuron.split,"/")
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_edgelist", folder = folder)
  gcsv = read.csv(gfile)
  gcsv
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_synapses <- function(savedir = TRUE, local = FALSE,
                            neuron.split = c("polypre_centrifugal_synapses",
                                             "polypre_centrifugal_distance")){
  neuron.split = match.arg(neuron.split)
  savedir = good_savedir(savedir = savedir,local = local)
  folder = paste0("hemibrain_neurons/",neuron.split,"/")
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_synapses", folder = folder)
  gcsv = read.csv(gfile)
  gcsv
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_connections <- function(savedir = TRUE, local = FALSE,
                               neuron.split = c("polypre_centrifugal_synapses",
                                                "polypre_centrifugal_distance")){
  neuron.split = match.arg(neuron.split)
  savedir = good_savedir(savedir = savedir,local = local)
  folder = paste0("hemibrain_neurons/",neuron.split,"/")
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_connections", folder = folder)
  gcsv = read.csv(gfile)
  gcsv
}


# hidden
find_gfile <- function(savedir,
                       file,
                       folder = "hemibrain_neurons/"){
  glist = list.files(paste0(savedir,folder),recursive = FALSE, full.names = TRUE)
  gfile = glist[grepl(file,glist)]
  gfile = sort(gfile)[1]
  gfile
}
