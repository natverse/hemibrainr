##############################################################################################
################################ Access Data on GoogleDrive ##################################
##############################################################################################

#' Read precomputed information from the hemibrain Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team Drive. This includes all synapses, neuron-neuron connections
#' and an edgelist for all hemibrain neurons, broken down by axon and dendrite assignments. NBLAST matrices for all neurons against all neurons in
#' rhe data set is also available, inluding ones broken down by neuron compartment.
#'
#' @param savedir directory in which neuron data has been deposited. If \code{TRUE} your default save directory is used, which is stored as: \code{options()$Gdrive_hemibrain_data}
#' @param local logical, whether to try to read locally saved neurons (by default at: \code{options()$hemibrain_data}) or neurons from Google Drive (\code{options()$Gdrive_hemibrain_data}).
#' @param neuron.split read saved neurons spit in which way? Folder names indicative of arguments passed to \code{\link{flow_centrality}}.
#' and \code{\link{hemibrain_flow_centrality}}
#' @param nblast the NBLAST you would like to retrieve, i.e. \code{"arbours"} gives you a normalised all by all NBLAST matrix of all branching arbour.
#' @param cable the type of cable we want to read. \code{"all"} indicates full neurons.
#' @param data the type of data to read, i.e. neurons, an nblast matrix or a \code{dotprops} object.
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
  gcsv = utils::read.csv(gfile)
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
  gcsv = utils::read.csv(gfile)
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
  gcsv = utils::read.csv(gfile)
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
  gcsv = utils::read.csv(gfile)
  gcsv
}

#' @export
#' @rdname hemibrain_googledrive
hemibrain_nblast <- function(savedir = TRUE,
                             local = FALSE,
                             nblast = c("all",
                                        "primary.neurites",
                                        "primary.dendrites",
                                        "axons",
                                        "dendrites",
                                        "spines",
                                        "tracts",
                                        "arbour",
                                        "simplified"),
                             neuron.split = c("polypre_centrifugal_synapses",
                                              "polypre_centrifugal_distance")){
  nblast = match.arg(nblast)
  neuron.split = match.arg(neuron.split)
  savedir = good_savedir(savedir = savedir,local = local)
  folder = if(nblast %in% c("all","spines")){
    paste0("hemibrain_nblast/")
  }else{
    paste0("hemibrain_nblast/nblast_",neuron.split,"/")
  }
  if(nblast=="all"){gfile = find_gfile(savedir = savedir, file = "hemibrain.aba.mean.compress.rda", folder = folder)}
  if(nblast=="spines"){gfile = find_gfile(savedir = savedir, file = "hemibrain.spine.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="primary.neurites"){gfile = find_gfile(savedir = savedir, file = "hemibrain.pnt.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="primary.dendrites"){gfile = find_gfile(savedir = savedir, file = "hemibrain.pd.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="axons"){gfile = find_gfile(savedir = savedir, file = "hemibrain.axon.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="dendrites"){gfile = find_gfile(savedir = savedir, file = "hemibrain.dendrite.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="arbour"){gfile = find_gfile(savedir = savedir, file = "hemibrain.arbour.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="tracts"){gfile = find_gfile(savedir = savedir, file = "hemibrain.tract.aba.mean.compressed.rda", folder = folder)}
  if(nblast=="simplified"){gfile = find_gfile(savedir = savedir, file = "hemibrain.simp.aba.mean.compressed.rda", folder = folder)}
  message("Loading NBLAST matrix from ", gfile)
  env <- new.env(parent = parent.frame())
  assign(nblast, get(load(gfile, env), envir = env), envir = env)
  return(env[[nblast]])
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

#' @export
#' @rdname hemibrain_googledrive
flycircuit_neurons <- function(savedir = TRUE,
                               local = FALSE,
                               folder = "hemibrain_neurons/light_level/flycircuit",
                               cable = c("all",
                                          "primary.neurites",
                                          "arbour"),
                               data = c("neuronlist",
                                        "nblast",
                                        "dps")){
  type = match.arg(type)
  data = match.arg(data)
  data = if(data!="neuronlist"){
    paste0("_",data)
  }else{
    ""
  }
  savedir = good_savedir(savedir = savedir,local = local)
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neuron_bodyids", folder = folder)
  if(type=="all"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_neurons%s.rds",data), folder = folder)}
  if(type=="primary.neurites"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_pnts.rds",data), folder = folder)}
  if(type=="all"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_arbour.rds",data), folder = folder)}
  readRDS(gfile)
}




