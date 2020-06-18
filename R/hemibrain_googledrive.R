##############################################################################################
################################ Access Data on GoogleDrive ##################################
##############################################################################################

#' Read precomputed information from the hemibrain Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team
#'   Drive. This includes all synapses, neuron-neuron connections and an
#'   edgelist for all hemibrain neurons, broken down by axon and dendrite
#'   assignments. NBLAST matrices for all neurons against all neurons in the
#'   data set are also available, including ones broken down by neuron
#'   compartment.
#'
#' @param savedir directory in which neuron data has been deposited. If
#'   \code{TRUE} your default save directory is used, which is stored as:
#'   \code{getOption("Gdrive_hemibrain_data")}
#' @param local logical, whether to try to read locally saved neurons (by
#'   default at: \code{getOption("hemibrain_data")}) or neurons from Google
#'   Drive (\code{getOption("Gdrive_hemibrain_data")}).
#' @param neuron.split read saved neurons split in which way? Folder names
#'   indicative of arguments passed to \code{\link{flow_centrality}}. and
#'   \code{\link{hemibrain_flow_centrality}}
#' @param nblast the NBLAST matrix you would like to retrieve, e.g.
#'   \code{"arbours"} gives you a normalised all by all NBLAST matrix of all
#'   branching arbour.
#' @param cable the type of cable we want to read. \code{"all"} indicates full
#'   neurons.
#' @param data the type of data to read, i.e. neurons, an NBLAST matrix or a
#'   \code{\link{dotprops}} object.
#' @param folder A subfolder containing the data object to read.
#'
#' @return a \code{data.frame} or character vector
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' #  All neuprint IDs for neurons that have a split precomputed
#' ids = hemibrain_neuron_bodyids()
#'
#' # Connectivity edgelist, broken down by axon/dendrite
#' elist = hemibrain_elist()
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}},
#'   \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{hemibrain_precomputed_splitpoints}},
#'   \code{\link{hemibrain_metrics}}
#' @name hemibrain_googledrive
#' @aliases hemibrain_neuron_bodyids
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
hemibrain_nblast <- function(nblast = c("all",
                                        "primary.neurites",
                                        "primary.dendrites",
                                        "axons",
                                        "dendrites",
                                        "spines",
                                        "tracts",
                                        "arbour",
                                        "simplified"),
                             savedir = TRUE,
                             local = FALSE,
                             neuron.split = c("polypre_centrifugal_synapses",
                                              "polypre_centrifugal_distance")){
  nblast = match.arg(nblast)
  neuron.split = match.arg(neuron.split)
  savedir = good_savedir(savedir = savedir,local = local)
  folder = if(nblast %in% c("all")){
    "hemibrain_nblast"
  }else{
    paste0("hemibrain_nblast/nblast_", neuron.split)
  }
  file = switch(nblast,
    all = "hemibrain.aba.mean.compress.rda",
    spines = "hemibrain.spine.aba.mean.compressed.rda",
    primary.neurites = "hemibrain.pnt.aba.mean.compressed.rda",
    primary.dendrites = "hemibrain.pd.aba.mean.compressed.rda",
    axons = "hemibrain.axon.aba.mean.compressed.rda",
    dendrites = "hemibrain.dendrite.aba.mean.compressed.rda",
    arbour = "hemibrain.arbour.aba.mean.compressed.rda",
    tracts = "hemibrain.tract.aba.mean.compressed.rda",
    simplified = "hemibrain.simp.aba.mean.compressed.rda",
    stop("Unrecognised value of nblast argument!")
  )
  gfile = find_gfile(savedir = savedir, file = file, folder = folder)
  message("Loading NBLAST matrix from ", gfile)
  env <- new.env(parent = parent.frame())
  assign(nblast, get(load(gfile, env), envir = env), envir = env)
  return(env[[nblast]])
}

# hidden
find_gfile <- function(savedir,
                       file,
                       folder = "hemibrain_neurons"){
  glist = list.files(file.path(savedir,folder), recursive = FALSE,
                     full.names = TRUE)
  gfile = glist[grepl(file, glist)]
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
                                        "dps"),
                               brainspace = c("JRCFIB2018F","JRCFIB2018Fraw","FCWB")){
  cable = match.arg(cable)
  data = match.arg(data)
  brainspace = match.arg(brainspace)
  brainspace = if(brainspace!="FCWB"){
    paste0("_",brainspace)
  }else{
    ""
  }
  if(data!="neuronlist"){
    data = paste0("_",data)
    brainspace = ""
  }else{
    data =  ""
  }
  savedir = good_savedir(savedir = savedir,local = local)
  if(cable=="all"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_neurons%s%s.rds",data,brainspace), folder = folder)}
  if(cable=="primary.neurites"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_pnts%s%s.rds",data,brainspace), folder = folder)}
  if(cable=="arbour"){gfile = find_gfile(savedir = savedir, file = sprintf("FlyCircuit_all_arbour%s.rds",data), folder = folder)}
  readRDS(gfile)
}

#' @export
#' @rdname hemibrain_googledrive
#' @param brainspace A template brain space for neurons loaded by
#'   \code{hemibrain_lm_lhns}. Defaults to \code{JRCFIB2018F}.
#' @importFrom utils installed.packages
hemibrain_lm_lhns <- function(savedir = TRUE,
                              local = FALSE,
                              folder = "hemibrain_neurons/light_level/lhns",
                              data = c("neuronlist",
                                       "nblast",
                                       "dps"),
                              cable = c("lhns","lhins","lines"),
                              brainspace = c("JRCFIB2018F","JRCFIB2018Fraw","FCWB")){
  brainspace = match.arg(brainspace)
  data = match.arg(data)
  cable = match.arg(cable)
  brainspace = if(brainspace=="FCWB") {
    # To be honest, this is a bit of a hack to avoid R CMD check
    # errors since lhns is not Suggested
    # I think lhns is too big and complicated a package
    if(isFALSE("lhns" %in% rownames(installed.packages())))
      stop("To use these data, please install the optional lhns package:\n",
           "remotes::install_github('jefferislab/lhns')")
    return(eval(str2lang("lhns::most.lhns")))
  } else {
    brainspace = paste0("_",brainspace)
  }
  if(data!="neuronlist"){
    data = paste0("_",data)
  }else{
    data =  ""
  }
  if(cable=="lines"){
    data = "_dps"
  }
  savedir = good_savedir(savedir = savedir,local = local)
  gfile = find_gfile(savedir = savedir, file = sprintf("most_%s%s%s.rds",cable,data,brainspace), folder = folder)
  readRDS(gfile)
}


