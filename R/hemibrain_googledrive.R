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
#' @param Gdrive name of the mounted Google drive. Can just be the name of a team
#' drive to be found in the standard location: \code{"/Volumes/GoogleDrive/Shared\ drives/"}.
#' @param path path to your 'hemibrainr Google drive', either a mounted remote or a local copy.
#' For example on a Mac, \href{https://support.google.com/drive/answer/7329379?authuser=2}{Google filestream} mounts it at: \code{/Volumes/GoogleDrive/Shared\ drives/}.
#' With \href{https://rclone.org/drive/}{rclone sheet} you can mount it wherever you like, by default hemibrainr expects it at:
#' \code{paste0(getwd(),"hemibrainr_rclone/")}
#'
#' @details In order to use this package to its fullest,
#' you need to get \code{hemibrainr} to read large amounts of data stored on a Google drive.
#' To do this, you have three options:
#'
#' 1. Mount your Google drive using \href{https://support.google.com/drive/answer/7329379?authuser=2}{Google filestream}
#' 2. Mount your Google drive using \href{https://rclone.org/drive/}{rclone}
#' 3. Download the Google drive and save locally
#'
#' The best option is to use Google filestream.
#' By default, this is what \code{hemibrainr} expects.
#' However, you need a \href{https://workspace.google.com/pricing.html}{Google Workspace}
#' account (formerly G-Suite), which is a paid-for service.
#' Without this, the best option is to use \href{https://rclone.org/drive/}{rclone}. For detailed
#' instructions on how to configure \code{rclone}, please see your article in this package:
#' Reading hemibrainr data from Google drive.
#' You can also download to an external hard drive and use that.
#'
#' We have two \href{https://support.google.com/a/users/answer/9310156?hl=en}{Google team drives}
#'  available for you to use, which contain similar data. One ("hemibrain")
#'  is for internal use by the \href{https://www.zoo.cam.ac.uk/research/groups/connectomics}{Drosophila Connectomics Group}.
#'  The other one (\code{"hemibrainr"}) is shared with those who would like access. Contact us by email to request access.
#'
#' @return Sets options \code{Gdrive_hemibrain_data}
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # The default drive is named 'hemibrainr'
#' # You could also set up your own drive with saved data.
#' hemibrainr_set_drive("hemibrainr")
#'
#' # Ah but what it you do not have the drive loaded
#' # With Google filestream say, as a network drive?
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
#' @seealso \code{\link{hemibrainr_googledrive_data}}
#' @name hemibrainr_set_drive
#' @export
hemibrainr_set_drive <- function(Gdrive = "hemibrainr",
                                 path = "/Volumes/GoogleDrive/Shared\ drives/"){
  if(!grepl("/",Gdrive)){
    Gdrive = sprintf("/Volumes/GoogleDrive/Shared\ drives/%s/", Gdrive)
  }
  options(Gdrive_hemibrain_data = paste0(path,Gdrive))
  message("Google drive path set to: ", options()$Gdrive_hemibrain_data)
  if(file.exists(options()$Gdrive_hemibrain_data)){
    message("Google drive found")
  }else{
    message("Google drive not found")
  }
}

#' @name hemibrainr_set_drive
#' @export
hemibrainr_team_drive <- function(){
  basename(options()$Gdrive_hemibrain_data)
}

#' Read precomputed hemibrain data from the hemibrainr Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team
#'   Drive. (see \code{\link{hemibrainr_set_drive}}) and (see \code{\link{hemibrainr_rclone}}).
#'   This includes body IDs for all hemibrain neurons ((\code{hemibrain_neuron_bodyids})),
#'   all synapses (\code{hemibrain_synapses}),
#'   neuron-neuron connections (\code{hemibrain_connections}) and an
#'   edgelist  (\code{hemibrain_elist}) for all hemibrain neurons, broken down by axon and dendrite
#'   assignments. Hemibrain related NBLASTs retrieved using \code{\link{hemibrain_nblast}}.
#'
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}),
#' but the user can specify an alternative path.
#' @param folder A subfolder on the hemibrain team drive or your local data folder
#'   containing the data object to read.
#' @param sql logical. Whether not to read the desired data from an \code{SQL} database. This can save you from having to load
#' a lot of information into memory, and it works well with \code{dplyr} pipes. If \code{FALSE} the relevant \code{.csv} is read from
#' the connected google drive (see \code{\link{hemibrainr_set_drive}}) and (see \code{\link{hemibrainr_rclone}}) and loaded into memory,
#' which can take some time.
#' @param ... if \code{sql=TRUE}, methods passed to \code{dplyr::tbl}.
#'
#' @return a \code{data.frame}. Depending on which synapse function was called, it can contain the columns:
#'
#' \itemize{
#'
#'   \item{"treenode_id"} { - the position of the node in the SWC-style table found at \code{neuron$d}, where the neuron is the skeleton for \code{bodyid}.}
#'
#'   \item{"connector_id"}{ - the unique ID for a pre/post synapse, as read from neuPrint. If this is not given, you are looking at a connection not a synapse.
#'   In this case \code{count} should be given, which shows the number of synapses in this connection.}
#'
#'   \item{"prepost"}{ - whether the given synapse is a pre-synape (0, output synapse) or postsynapse (1, input synapse). Alternatively, if a connection is given,
#'   whether this connection is presynaptic to \code{bodyid} (0, \code{bodyid} is target) or postsynaptic (1, \code{bodyid} is source).}
#'
#'   \item{"x"}{ - x coordinate for the root point.}
#'
#'   \item{"y"}{ - y coordinate for the root point.}
#'
#'   \item{"z"}{ - z coordinate for the root point.}
#'
#'   \item{"confidence"}{ - FlyEM's confidence level. The lower the score, the more likely this synapse is an artefact.}
#'
#'   \item{"bodyid"}{ - The neuPrint neuron/body related to the synapse/connection given in each row.}
#'
#'   \item{"partner"}{ - The neuron connecting to \code{bodyid} by the givne synapse/connection.}
#'
#'   \item{"pre"}{ - The body ID for the presynaptic (source) neuron.}
#'
#'   \item{"partner"}{ - The body ID for the presynaptic (target) neuron.}
#'
#'   \item{"Label"}{ - The compartment of the \code{bodyid} neuron on which the synapse is placed / which receives/makes the given connection.
#'   See \code{?standardise}.}
#'
#'   \item{"partner.Label"}{ - The compartment of the \code{partner} neuron on which the synapse is placed / which receives/makes the given connection.}
#'
#'   \item{"count"}{ - The number of synapses that make the given connection. Sometimes referred to as 'weight'.}
#'
#'   \item{"norm"}{ - The normalised synapse weight. \code{count} is divided by the total number of inputs that the
#'   target neuron's (\code{post}) compartment (\code{Label}) has. I.e. this normalisation is by total inputs onto a dendrite or axon, not the whole neuron.}
#'
#'}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # All neuprint IDs for neurons that have a split precomputed
#' ids = hemibrain_neuron_bodyids()
#'
#' # For these body IDs, all synapses:
#' syns = hemibrain_synapses()
#'
#' # For these body IDs, input and output connections:
#' conns = hemibrain_connections()
#'
#' # Connectivity edgelist, broken down by axon/dendrite
#' elist = hemibrain_elist()
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}},
#'   \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{hemibrainr_googledrive_data}},
#'   \code{\link{hemibrain_metrics}}
#' @name hemibrainr_googledrive_data
#' @aliases hemibrain_neuron_bodyids
#' @export
hemibrain_neuron_bodyids <- function(local = FALSE, folder = "hemibrain_neurons/", sql = FALSE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "hemibrain_all_neuron_bodyids", sql.db = "hemibrainr_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neuron_bodyids", folder = folder)
    gcsv = as.data.frame(readr::read_csv(gfile))
    as.character(gcsv$x)
  }
}

#' @rdname hemibrainr_googledrive_data
#' @export
hemibrain_elist <- function(local = FALSE, folder = "hemibrain_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    gcsv = find_gsql(savedir = savedir, tab = "hemibrain_all_neurons_edgelist_polypre_centrifugal_synapses", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_edgelist", folder = folder)
    gcsv = as.data.frame(readr::read_csv(gfile))
  }
  gcsv
}

#' @rdname hemibrainr_googledrive_data
#' @export
hemibrain_synapses <- function(local = FALSE, folder = "hemibrain_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    gcsv = find_gsql(savedir = savedir, tab = "hemibrain_all_neurons_synapses_polypre_centrifugal_synapses", sql.db = "hemibrainr_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_synapses", folder = folder)
    gcsv = as.data.frame(readr::read_csv(gfile))
  }
}

#' @rdname hemibrainr_googledrive_data
#' @export
hemibrain_connections <- function(local = FALSE, folder = "hemibrain_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    gcsv = find_gsql(savedir = savedir, tab = "hemibrain_all_neurons_connections_polypre_centrifugal_synapses", sql.db = "hemibrainr_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "hemibrain_all_neurons_connections", folder = folder)
    gcsv = as.data.frame(readr::read_csv(gfile))
  }
}

#' Read precomputed NBLASTs from the hemibrainr Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team
#'   Drive. This includes body IDs for all hemibrain neurons ((\code{hemibrain_neuron_bodyids})),
#'   all synapses (\code{hemibrain_synapses}),
#'   neuron-neuron connections (\code{hemibrain_connections}) and an
#'   edgelist  (\code{hemibrain_elist}) for all hemibrain neurons, broken down by axon and dendrite
#'   assignments. NBLAST matrices for all neurons against all neurons in the
#'   data set are also available, including ones broken down by neuron
#'   compartment.
#'
#' @param nblast the NBLAST matrix you would like to retrieve, e.g.
#'   \code{"arbours"} gives you a normalised all by all NBLAST matrix of all
#'   branching arbour.
#'   \code{\link{dotprops}} object.
#' @inheritParams hemibrainr_googledrive_data
#'
#' @details NBLASTs were made in \code{JRCFIB2018F} space. Hemibrain neurons first had their 'twigs' pruned:
#'
#' \code{fib.twigs5 = nlapply(all.neurons.flow, prune_twigs, twig_length=5000, .parallel = TRUE, OmitFailures = TRUE)}.
#'
#' Where \code{all.neurons.flow} can be called with \code{hemibrain_neurons}. Neurons were moved from \code{JRC2018Fraw}
#' to \code{JRC2018F} by:
#'
#' \code{all.neurons.flow.microns = hemibrainr:::scale_neurons.neuronlist(fib.twigs5, .parallel = TRUE, OmitFailures = TRUE)}
#'
#' For 'compartment' NBLASTs neurons' axons, dendrites, primary neurite
#' tracts, etc., were extracted. E.g.
#'
#' \code{all.neurons.pnt = primary_neurite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}
#' \code{all.neurons.arbour = arbour_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}
#' \code{all.neurons.tract = tract_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}
#'
#' For the 'simplified' hemibrain neuron NBLAST, the function \code{nat::simplify_neuron} was used as so:
#'
#' \code{all.neurons.simp = nat::nlapply(X = all.neurons.flow.microns, FUN = nat::simplify_neuron, n = 1, invert = FALSE, .parallel = TRUE, OmitFailures = TRUE)}
#'
#' Neurons, simplified neurons and compartments were turned into vector cloud using \code{nat::dotprops} as so:
#'
#' \code{hemibrain.microns.dps =dotprops(all.neurons.flow.microns, k=5, resample=1, .parallel=T, OmitFailures = T)}
#'
#' NBLASTs were run as so:
#'
#' \code{hemibrain.aba.mean=nat.nblast::nblast_allbyall(hemibrain.twigs5.dps,
#' .progress='text',.parallel=TRUE, normalisation='mean')}
#'
#' NBLASTs were also run against flywire and flycircuit neurons.
#' Flywire neurons were pulled from the \href{https://ngl.flywire.ai/?local_id=1191adb1d6f85f5f8be0fedc675460f6}{flywire project} as meshes and skeletonised
#' using \code{fafbseg::skeletor}.
#'
#' The code for this pipeline can be found as a \href{https://github.com/flyconnectome/fafbpipeline}{flyconnectome GitHub repository}
#'
#' @return a \code{matrix} with named rows and columns describing an NBLAST result.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # A normalised neuron-neuron NBLASST for all hemibrain neurons
#' hemibrain.nblast = hemibrain_nblast(nblast = "hemibrain")
#'
#' # A normalised neuron-neuron NBLASST for all hemibrain neurons's axons
#' hemibrain.axon.nblast = hemibrain_nblast(nblast = "hemibrain-axons")
#'
#' # And for simplified hemibrain neurons
#' hemibrain.simp.nblast = hemibrain_nblast(nblast = "hemibrain-simplified")
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}},
#'   \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{hemibrain_metrics}},
#'   \code{\link{hemibrainr_googledrive_data}}
#' @name hemibrain_nblast
#' @export
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

# hidden
find_gsql <- function(savedir,
                      sql.db = "hemibrainr_data.sqlite",
                      tab = NULL,
                      folder = "hemibrain_neurons",
                      ...){
  gfile = find_gfile(savedir = savedir, file = sql.db, folder = folder)
  sql.db.load  = dplyr::src_sqlite(gfile, create = FALSE)
  if(!is.null(tab)){
    dplyr::tbl(sql.db.load, tab, ...)
  }else{
    sql.db.load
  }
}

#' Read flycircuit neurons from hemibrainr Google Drive
#'
#' @description Read \href{http://www.flycircuit.tw/}{flycircuit} neurons from the \code{hemibrainr} google drive.
#' The google drive must be mounted with Google Filestream or rclone.
#' See \code{\link{hemibrainr_set_drive}}.
#'
#' @param cable the type of cable we want to read. \code{"all"} indicates full neurons.
#' @param data the type of data to read, i.e. neurons, an NBLAST matrix or a
#'   \code{\link{dotprops}} object.
#' @param brainspace A template brain space for neurons loaded by
#'   \code{lm_lhns}. Defaults to \code{JRCFIB2018F}.
#' @inheritParams hemibrainr_googledrive_data
#'
#' @return a \code{nat::neuronlistfh} object for 'light-level' neurons in the given brainspace.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # All flycircuit neurons
#' fw = flycircuit_neurons()
#'
#' # All LHNs from Frechter et al. 2019, eLife
#' lhns = lm_lhns()
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}},
#'   \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{hemibrain_metrics}},
#'   \code{\link{hemibrainr_googledrive_data}}
#' @name flycircuit_neurons
#' @export
flycircuit_neurons <- function(local = FALSE,
                               folder = "light_level/flycircuit",
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
#' @rdname flycircuit_neurons
#' @importFrom utils installed.packages
lm_lhns <- function(local = FALSE,
                              folder = "light_level/lhns",
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


