##############################################################################################
################################ Access Data on GoogleDrive ##################################
##############################################################################################

#' Set data source for precomputed information
#'
#' @description One can read precomputed data available on the hemibrain Google Team
#'   Drive, and other supported team drives. This includes all synapses, neuron-neuron connections and an
#'   edgelist for all hemibrain neurons, broken down by axon and dendrite
#'   assignments. NBLAST matrices for all neurons against all neurons in the
#'   data set are also available, including ones broken down by neuron
#'   compartment.
#'
#' @param Gdrive name of the mounted google drive. Can just be the name of a team
#' drive to be found in the standard location: \code{"/Volumes/GoogleDrive/Shared\ drives/"}.
#' @param path path to your 'hemibrainr google drive', either a mounted remote or a local copy.
#' For example on a Mac, \href{https://support.google.com/drive/answer/7329379?authuser=2}{Google filestream} mounts it at: \code{/Volumes/GoogleDrive/Shared\ drives/}.
#' With \href{https://rclone.org/drive/}{rclone sheet} you can mount it wherever you like, by default hemibrainr expects it at:
#' \code{paste0(getwd(),"hemibrainr_rclone/")}
#'
#' @details In order to use this package to its fullest,
#' you need to get \code{hemibrainr} to read large amounts of data stored on a google drive.
#' To do this, you have three options:
#'
#' 1. Mount your google drive usin \href{https://support.google.com/drive/answer/7329379?authuser=2}{Google filestream}
#' 2. Mount your google drive using \href{https://rclone.org/drive/}{rclone}
#' 3. Download the google drive and save locally
#'
#' The best option is to use google filestream.
#' By default, this is what \code{hemibrainr} expects.
#' However, you need a \href{https://workspace.google.com/pricing.html}{Google Workspace}
#' account (formerly G-Suite), which is a paid-for service.
#' Without this, the best option is to use href{https://rclone.org/drive/}{rclone}. For detailed
#' instructions on how to confiure rclone, please see your article in this package:
#' Reading hemibrainr data from google drive.
#' You can also download to an external hard drive and use that.
#'
#' We have two \href{https://support.google.com/a/users/answer/9310156?hl=en}{Google team drives}
#'  available for you to use, which contain similar data. One ("hemibrain")
#'  is for internal use by the \href{https://www.zoo.cam.ac.uk/research/groups/connectomics}{Drosphila Connctomcis Group}.
#'  The other one (\code{"hemibrainr"}) is shared with those who would like access. Contact us by email to request access.
#'
#' @return Sets options \code{Gdrive_hemibrain_data}
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # The default drive is named 'hemibrain'
#' # The Wilson lab at HMS uses a drive called 'hemibrainr'.
#' # You could also set up your own drive with saved data.
#' hemibrainr_set_drive("hemibrainr")
#'
#' # Ah but what it you do not have the dirve loaded
#' # With google filestream say, as a network drive?
#' # In this case, you can use the free-service rclone instead
#' # However, you will need to install and configure it:https://rclone.org/drive/
#' hemibrainr_rclone(Gdrive = "hemibrainr")
#' ### Mounts in your working directory
#'
#' # All neuprint IDs for neurons that have a split precomputed
#' ids = hemibrain_neuron_bodyids()
#'
#' # Connectivity edgelist, broken down by axon/dendrite
#' elist = hemibrain_elist()
#'
#' # Unmount rclone drive
#' hemibrainr_rclone_unmount()
#'
#' }}
#' @seealso \code{\link{hemibrain_googledrive}}
#' @name hemibrainr_set_drive
#' @export
hemibrainr_set_drive <- function(Gdrive = "hemibrain",
                                 path = "/Volumes/GoogleDrive/Shared\ drives/"){
  if(!grepl("/",Gdrive)){
    Gdrive = sprintf("/Volumes/GoogleDrive/Shared\ drives/%s/", Gdrive)
  }
  options(Gdrive_hemibrain_data = paste0(path,Gdrive))
  message("Google drive path set to: ", options()$Gdrive_hemibrain_data)
  if(file.exists(options()$Gdrive_hemibrain_data)){
    message("Google drive found")
  }else{
    message("Google not found")
  }
}

#' @name hemibrainr_set_drive
#' @export
hemibrainr_team_drive <- function(){
  basename(options()$Gdrive_hemibrain_data)
}

#' Read precomputed information from the hemibrain Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team
#'   Drive. This includes all synapses, neuron-neuron connections and an
#'   edgelist for all hemibrain neurons, broken down by axon and dendrite
#'   assignments. NBLAST matrices for all neurons against all neurons in the
#'   data set are also available, including ones broken down by neuron
#'   compartment.
#'
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}),
#' but the user can specit and alternative path.
#' @param nblast the NBLAST matrix you would like to retrieve, e.g.
#'   \code{"arbours"} gives you a normalised all by all NBLAST matrix of all
#'   branching arbour.
#' @param cable the type of cable we want to read. \code{"all"} indicates full
#'   neurons.
#' @param data the type of data to read, i.e. neurons, an NBLAST matrix or a
#'   \code{\link{dotprops}} object.
#' @param folder A subfolder on the Hemibrain team drive or your local data foldercontaining the data object to read.
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
hemibrain_neuron_bodyids <- function(local = FALSE, folder = "hemibrain_neurons/"){
  savedir = good_savedir(local = local)
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neuron_bodyids", folder = folder)
  gcsv = utils::read.csv(gfile)
  as.character(gcsv$x)
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_elist <- function(local = FALSE, folder = "hemibrain_neurons/"){
  savedir = good_savedir(local = local)
  folder = "hemibrain_neurons/"
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_edgelist", folder = folder)
  gcsv = utils::read.csv(gfile)
  gcsv
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_synapses <- function(local = FALSE, folder = "hemibrain_neurons/"){
  savedir = good_savedir(local = local)
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_synapses", folder = folder)
  gcsv = utils::read.csv(gfile)
  gcsv
}

#' @rdname hemibrain_googledrive
#' @export
hemibrain_connections <- function(local = FALSE, folder = "hemibrain_neurons/"){
  savedir = good_savedir(local = local)
  gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_connections", folder = folder)
  gcsv = utils::read.csv(gfile)
  gcsv
}

#' @export
#' @rdname hemibrain_googledrive
hemibrain_nblast <- function(nblast = c("hemibrain",
                                        "flywire",
                                        "hemibrain-flywire",
                                        "hemibrain-fafb14",
                                        "hemibrain-flycircuit",
                                        "flywire-mirror",
                                        "hemibrain-primary.neurites",
                                        "hemibrain-primary.dendrites",
                                        "hemibrain-axons",
                                        "hemibrain-dendrites",
                                        "hemibrain-spines",
                                        "hemibrain-tracts",
                                        "hemibrain-arbour",
                                        "hemibrain-simplified"),
                             local = FALSE){
  nblast = match.arg(nblast)
  savedir = good_savedir(local = local)
  folder = if(nblast %in% c("hemibrain",
                            "hemibrain-flywire",
                            "hemibrain-flycircuit",
                            "hemibrain-fafb14",
                            "flywire-mirror")){
    "hemibrain_nblast"
  }else{
    "hemibrain_nblast/nblast_"
  }
  file = switch(nblast,
    hemibrain = "hemibrain.aba.mean.compress.rda",
    `hemibrain-spines` = "hemibrain.spine.aba.mean.compressed.rda",
    `hemibrain-primary.neurites` = "hemibrain.pnt.aba.mean.compressed.rda",
    `hemibrain-primary.dendrites` = "hemibrain.pd.aba.mean.compressed.rda",
    `hemibrain-axons` = "hemibrain.axon.aba.mean.compressed.rda",
    `hemibrain-dendrites` = "hemibrain.dendrite.aba.mean.compressed.rda",
    `hemibrain-arbour` = "hemibrain.arbour.aba.mean.compressed.rda",
    `hemibrain-tracts` = "hemibrain.tract.aba.mean.compressed.rda",
    `hemibrain-simplified` = "hemibrain.simp.aba.mean.compressed.rda",
    `hemibrain-flywire` = "hemibrain.flywire.mean.compressed.rda",
    `hemibrain-flycircuit` = "hemibrain.flycircuit.mean.compressed.rda",
    `hemibrain-fafb14` = "hemibrain.fafb14.mean.compressed.rda",
    `flywire-mirror` = "flywire.mirror.mean.compressed.rda",
    stop("Unrecognised value of nblast argument!")
  )
  if(nblast=="hemibrain-fafb14" && file.exists("/Volumes/GoogleDrive/Shared drives/flyconnectome/fafbpipeline/fib.fafb.crossnblast.twigs5.mean.compress.rda")){
    gfile = "/Volumes/GoogleDrive/Shared drives/flyconnectome/fafbpipeline/fib.fafb.crossnblast.twigs5.mean.compress.rda"
  }else{
    gfile = find_gfile(savedir = savedir, file = file, folder = folder)
  }
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
flycircuit_neurons <- function(local = FALSE,
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
  savedir = good_savedir(local = local)
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
hemibrain_lm_lhns <- function(local = FALSE,
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
  savedir = good_savedir(local = local)
  gfile = find_gfile(savedir = savedir, file = sprintf("most_%s%s%s.rds",cable,data,brainspace), folder = folder)
  readRDS(gfile)
}


